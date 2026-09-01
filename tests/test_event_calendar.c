#include "EDMD.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

extern int N;
extern int actualPaulList;
extern int paulListN;
extern double dtPaul;
extern double paulTime;
extern node *root;
extern node **eventPaul;

static int failures;
static node sentinel;
static node *buckets[65];

static void fail(const char *message)
{
	fprintf(stderr, "%s\n", message);
	failures++;
}

static void expect_true(int condition, const char *message)
{
	if (!condition)
		fail(message);
}

static void expect_close(double actual, double expected, const char *message)
{
	const double scale = fmax(1.0, fmax(fabs(actual), fabs(expected)));
	if (fabs(actual - expected) > 1e-13*scale)
		fail(message);
}

static node blank_event(double eventTime)
{
	node result;
	memset(&result, 0, sizeof(result));
	result.t = eventTime;
	result.q = -1;
	return result;
}

/* Each test drains the preceding queue before calling this helper. */
static void reset_calendar(int bucketCount, double bucketWidth,
	double origin, int currentBucket)
{
	N = bucketCount;
	paulListN = bucketCount;
	dtPaul = bucketWidth;
	paulTime = origin;
	actualPaulList = currentBucket;
	memset(buckets, 0, sizeof(buckets));
	eventPaul = buckets;
	memset(&sentinel, 0, sizeof(sentinel));
	sentinel.t = never + 1.0;
	sentinel.q = -1;
	root = &sentinel;
}

static void expect_detached(const node *eventNode, const char *message)
{
	expect_true(eventNode->q == -1 && eventNode->lft == NULL &&
		eventNode->rgt == NULL && eventNode->top == NULL, message);
}

static void test_invalid_admission_and_unique_requeue(void)
{
	reset_calendar(4, 1.0, 0.0, 0);
	double invalidTimes[] = {NAN, INFINITY, -INFINITY, never, 2.0*never};
	for (size_t index = 0;
		index < sizeof(invalidTimes)/sizeof(invalidTimes[0]); index++){
		node invalid = blank_event(invalidTimes[index]);
		addEventToQueue(&invalid);
		expect_detached(&invalid, "invalid event retained queue links");
	}
	expect_true(findNextEvent() == NULL,
		"invalid events made the calendar non-empty");

	node replacement = blank_event(0.5);
	addEventToQueue(&replacement);
	replacement.t = 0.75;
	addEventToQueue(&replacement);
	expect_true(findNextEvent() == &replacement,
		"requeued fixed-storage node was not the unique minimum");
	removeEventFromQueue(&replacement);
	expect_detached(&replacement, "removed tree event retained queue links");
	expect_true(findNextEvent() == NULL,
		"duplicate fixed-storage node remained after removal");
}

static void test_horizon_rounding_does_not_alias_current_bucket(void)
{
	const int bucketCount = 39;
	const double bucketWidth = 5.0/bucketCount;
	reset_calendar(bucketCount, bucketWidth, 0.0, 0);
	node edge = blank_event(nextafter(bucketWidth*bucketCount, -INFINITY));
	addEventToQueue(&edge);
	expect_true(edge.q > 0 && edge.q < bucketCount,
		"near-horizon event aliased the current or overflow bucket");
	removeEventFromQueue(&edge);
	expect_detached(&edge, "near-horizon removal corrupted node links");
	expect_true(findNextEvent() == NULL,
		"near-horizon event remained scheduled after removal");
}

static void test_empty_bucket_and_sparse_epoch_jumps(void)
{
	reset_calendar(8, 0.25, 0.0, 0);
	node laterBucket = blank_event(1.55);
	addEventToQueue(&laterBucket);
	expect_true(findNextEvent() == &laterBucket,
		"calendar did not find the next occupied regular bucket");
	expect_close(paulTime, 1.5,
		"calendar did not jump directly across empty regular buckets");
	removeEventFromQueue(&laterBucket);
	expect_true(findNextEvent() == NULL,
		"regular-bucket jump left a stale event");

	reset_calendar(4, 1.0, 0.0, 0);
	node farFuture = blank_event(1e20);
	addEventToQueue(&farFuture);
	expect_true(farFuture.q == paulListN,
		"far-future event was not placed in overflow");
	expect_true(findNextEvent() == &farFuture,
		"sparse overflow event was not recovered");
	expect_true(paulTime == farFuture.t,
		"sparse calendar walked epochs instead of rebasing once");
	removeEventFromQueue(&farFuture);
	expect_true(findNextEvent() == NULL,
		"sparse overflow event remained after removal");
}

static void test_wrap_redistributes_overflow_before_later_buckets(void)
{
	reset_calendar(4, 1.0, 0.0, 0);
	node overflow = blank_event(4.1);
	node marker = blank_event(2.1);
	addEventToQueue(&overflow);
	addEventToQueue(&marker);
	expect_true(findNextEvent() == &marker,
		"calendar did not reach the marker before the first wrap");
	expect_close(paulTime, 2.0,
		"pre-wrap marker bucket started at the wrong time");
	removeEventFromQueue(&marker);

	/* The overflow node has now aged to 2.1 time units away.  A newly queued
	 * regular node wraps to bucket zero, but the overflow node is still earlier
	 * and must be redistributed before that bucket is selected. */
	node afterWrap = blank_event(4.5);
	addEventToQueue(&afterWrap);
	expect_true(findNextEvent() == &overflow,
		"aged overflow did not preempt a later wrapped regular event");
	expect_close(paulTime, 4.0,
		"overflow redistribution occurred at the wrong wrap boundary");
	removeEventFromQueue(&overflow);
	expect_true(findNextEvent() == &afterWrap,
		"wrapped regular event was lost after overflow removal");
	removeEventFromQueue(&afterWrap);
	expect_true(findNextEvent() == NULL,
		"cross-wrap test did not drain the queue");
}

int main(void)
{
	test_invalid_admission_and_unique_requeue();
	test_horizon_rounding_does_not_alias_current_bucket();
	test_empty_bucket_and_sparse_epoch_jumps();
	test_wrap_redistributes_overflow_before_later_buckets();

	if (failures != 0){
		fprintf(stderr, "%d event-calendar regression(s) failed\n", failures);
		return 1;
	}
	printf("event-calendar regressions passed\n");
	return 0;
}
