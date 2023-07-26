# Graphical Event Driven Molecular Dynamics

Usually, to deal with collisions, one integrates the equation of motion using a small timestep. A notoriously hard task is to deal with object interpenetration. Tricks must be used to supress these overlaps and the resulting dynamics is **always** nonphysical. Artifacts due to these corrections often prevent high densities to be reached. Here, (almost) everything is solved exactly by computing the time before the next collision happening in the system and then advancing the simulation time to this next collision + processing it.
This type of algorithm is extremely fast but you pay the price of having to deal with an asynchronized evolution of events. Since every action is driven by the scheduling of an "event" (every collision in the future, every cell partitioning, every display of the screen in the game loop becomes an event), one has to implement a fast data structure able to quickly sort these events, retrieve them, delete them, reschedule them. This is achieved by a BST and a buffer queue.

At full speed, and for ~2000 particles the algorithm can perform around $10^6$ collisions/sec.

## Compiling

For now, the project can only be compiled on linux but it should be pretty straighforward to modify it to compile on Windows. The graphical part probably only works with x11 but it's a matter of finding the good makefile :).

Use ```make``` for the algorithm with the graphical interface and ```make G=0``` for the algorithm without a graphical interface.

##### Dependencies:
  For the graphical part: [raylib](https://github.com/raysan5/raylib/wiki/Working-on-GNU-Linux)

## Running

Compile and execute the ```./a.out```. You can pass option in the command line such as ```./a.out -N 1000 --phi 0.8```. This is particularly useful for the no-graphics mode of the program but almost useless when the graphical interface is used.
Without graphical interface the program ouputs files readable by [ovito](https://www.ovito.org/about/).
## Structure of the program

#### HUGE DISCLAIMER: I AM NOT --BY ANY STANDARD-- A PROGRAMMER. DO NOT USE MY PROGRAM AS A STANDARD FOR THIS KIND OF SIMULATION/GRAPHICAL DISPLAY
### Technical part
Everything is pure C no library.

- Particles are stored in cell-lists by an array based doubly linked list in order to efficiently compute the time before the next collision.
- Events (particles crossing cell-lists boundaries, particles colliding, thermostating, displaying, taking user inputs) are stored in a Binary Search Tree which receive events from [a partially sorted buffer/priority queue](https://arxiv.org/abs/physics/0606226).

### Graphical part
This program was first designed without a graphical interface in mind and uses a lot of global variable. The graphical interface **must not** interfere in any way to the structure of the initial program. Which might partially explains some poor design choice.
If it has to be redone, I would change a lot of things and encapsulate more the code. I'm still reworking the graphical part to alleviate the use of global variables.

## Sources

If you want to understand the structure of an EDMD algorithm: [Github repo EDMD](https://github.com/FSmallenburg/EDMD) with [research article](https://arxiv.org/abs/2201.01100).
