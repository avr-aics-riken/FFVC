#ifndef tt_time_h
#define tt_time_h

#if defined(WIN32)
// for Windows
#include <time.h>
class TtTime {
public:
	TtTime() {}

	void start() { clock_start = clock(); }
	void end()   { clock_end   = clock(); }

	// return elapsed time [seconds]
	float getElapsedSec() const {
		return float(clock_end - clock_start) / CLOCKS_PER_SEC;
	}

	// return elapsed time [milliseconds]
	float getElapsedMSec() const {
		return 1000.f * float(clock_end - clock_start) / CLOCKS_PER_SEC;
	}

private:
	clock_t clock_start, clock_end;
};
#else
// for UNIX
#include <sys/time.h>
class TtTime {
public:
	TtTime() {}

	void start() { gettimeofday(&tv_start, 0); }
	void end()   { gettimeofday(&tv_end,   0); }

	// return elapsed time [seconds]
	float getElapsedSec() const {
		return (tv_end.tv_sec - tv_start.tv_sec)
			+ .000001 * (tv_end.tv_usec - tv_start.tv_usec);
	}

	// return elapsed time [milliseconds]
	float getElapsedMSec() const {
		return 1000. * (tv_end.tv_sec - tv_start.tv_sec)
			+ .001 * (tv_end.tv_usec - tv_start.tv_usec);
	}

private:
	struct timeval tv_start, tv_end;
};
#endif

#endif  // tt_time_h

