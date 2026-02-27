import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from matplotlib.patches import Polygon as MplPolygon, Circle as MplCircle


# ----------------------------
# Rigid bodies
# ----------------------------
@dataclass
class Circle:
    pos: np.ndarray      # shape (2,)
    vel: np.ndarray      # shape (2,)
    radius: float

    def center(self, t: float) -> np.ndarray:
        return self.pos + self.vel * t


@dataclass
class RigidPolygon:
    verts_body: np.ndarray  # shape (N,2), vertices in the body frame (about origin)
    pos: np.ndarray         # shape (2,)
    vel: np.ndarray         # shape (2,)
    theta: float            # radians
    omega: float            # rad/s

    def verts_world(self, t: float) -> np.ndarray:
        th = self.theta + self.omega * t
        c, s = np.cos(th), np.sin(th)
        R = np.array([[c, -s],
                      [s,  c]])
        # row vectors: v_world = v_body @ R^T
        return (self.pos + self.vel * t) + (self.verts_body @ R.T)


# ----------------------------
# Geometry helpers
# ----------------------------
def point_segment_distance(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    """Euclidean distance from point p to segment [a,b]."""
    ab = b - a
    denom = float(np.dot(ab, ab))
    if denom == 0.0:
        return float(np.linalg.norm(p - a))
    t = float(np.dot(p - a, ab) / denom)
    t = float(np.clip(t, 0.0, 1.0))
    closest = a + t * ab
    return float(np.linalg.norm(p - closest))


def point_in_polygon(p: np.ndarray, poly: np.ndarray) -> bool:
    """
    Ray casting test for point in polygon (works for simple polygons, convex or concave).
    Note: boundary cases can be ambiguous; fine for a toy prototype.
    """
    x, y = float(p[0]), float(p[1])
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = float(poly[i][0]), float(poly[i][1])
        x2, y2 = float(poly[(i + 1) % n][0]), float(poly[(i + 1) % n][1])
        if (y1 > y) != (y2 > y):
            xinters = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < xinters:
                inside = not inside
    return inside


def signed_distance_point_polygon(p: np.ndarray, poly: np.ndarray) -> float:
    """
    Signed distance from point to polygon boundary:
      > 0 outside, < 0 inside.
    Uses min distance to edges + point-in-poly for the sign.
    """
    dmin = float("inf")
    n = len(poly)
    for i in range(n):
        a = poly[i]
        b = poly[(i + 1) % n]
        dmin = min(dmin, point_segment_distance(p, a, b))
    return -dmin if point_in_polygon(p, poly) else dmin


# ----------------------------
# Continuous collision detection (prototype)
# ----------------------------
def collision_time_circle_polygon(
    circle: Circle,
    poly: RigidPolygon,
    t_max: float = 10.0,
    dt: float = 1e-3,
    tol: float = 1e-10,
    max_bisect_iters: int = 80,
):
    """
    Find earliest collision time t in [0, t_max] where:
        signed_distance(circle_center(t), polygon(t)) - radius == 0
    using:
      1) time stepping to bracket a sign change
      2) bisection refinement

    Returns:
      (t_hit, info_dict) or (None, info_dict)
    """

    def g(t: float) -> float:
        p = circle.center(t)
        verts = poly.verts_world(t)
        return signed_distance_point_polygon(p, verts) - circle.radius

    g0 = g(0.0)
    if g0 <= 0.0:
        return 0.0, {"status": "initially_overlapping_or_touching", "g0": g0}

    steps = int(np.ceil(t_max / dt))
    t_prev, g_prev = 0.0, g0

    for k in range(1, steps + 1):
        t = min(k * dt, t_max)
        g_now = g(t)

        # Collision happened when crossing from positive to <= 0
        if g_now <= 0.0:
            a, b = t_prev, t  # bracket
            # Bisection refinement
            for _ in range(max_bisect_iters):
                m = 0.5 * (a + b)
                gm = g(m)
                if gm > 0.0:
                    a = m
                else:
                    b = m
                if (b - a) <= tol:
                    break

            return b, {
                "status": "hit",
                "bracket": (t_prev, t),
                "dt": dt,
                "tol": tol,
            }

        t_prev, g_prev = t, g_now
        if t >= t_max:
            break

    return None, {"status": "no_hit_in_interval", "t_max": t_max, "dt": dt}


# ----------------------------
# Visualization
# ----------------------------
def plot_snapshot(ax, circle: Circle, poly: RigidPolygon, t: float, title: str):
    verts = poly.verts_world(t)
    c = circle.center(t)

    ax.add_patch(MplPolygon(verts, fill=False))
    ax.add_patch(MplCircle(c, circle.radius, fill=False))

    pts = np.vstack([verts, c[None, :]])
    pad = 1.5 * circle.radius + 0.5
    mn = pts.min(axis=0) - pad
    mx = pts.max(axis=0) + pad
    ax.set_xlim(float(mn[0]), float(mx[0]))
    ax.set_ylim(float(mn[1]), float(mx[1]))
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title)
    ax.grid(True)


# ----------------------------
# Demo / sanity check
# ----------------------------
if __name__ == "__main__":
    # Body-frame polygon (convex-ish pentagon)
    verts_body = np.array([
        [-0.7, -0.3],
        [ 0.2, -0.3],
        [ 0.2,  0.0],
        [ -0.72,  0.7],
        [-0.6,  0.4],
    ], dtype=float)

    circle = Circle(
        pos=np.array([-3.0, 3]),
        vel=np.array([1, -1]),
        radius=0.1
    )

    poly = RigidPolygon(
        verts_body=verts_body,
        pos=np.array([0.0, 0.0]),
        vel=np.array([0.0, -0.05]),
        theta=0.5,
        omega=10.9
    )
    import time
    
    time_start = time.time()
    t_hit, info = collision_time_circle_polygon(circle, poly, t_max=50.0, dt=2e-1, tol=1e-10)
    time_end = time.time()
    print(f"Collision time computation took {time_end - time_start:.4f} seconds")
    print("collision_time:", t_hit)
    print("info:", info)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    plot_snapshot(axes[0], circle, poly, 0.0, "t = 0 (initial)")

    if t_hit is None:
        axes[1].set_title("No collision in interval")
        axes[1].axis("off")
    else:
        plot_snapshot(axes[1], circle, poly, t_hit, f"t = {t_hit:.6f} (collision)")

    fig.tight_layout()
    plt.show()