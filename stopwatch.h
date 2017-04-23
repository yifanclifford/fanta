#include <sys/time.h>
#include <stdlib.h>

class Stopwatch
{
public:
	float total_time;
	timeval t1, t2;
	Stopwatch(): total_time(0) {}
	float get_time() { return total_time; }	// 0 before the first stop

	void start() {
		gettimeofday(&t1, NULL);
	}

	void resume() {
		gettimeofday(&t1, NULL);
	}

	void pause() {
		gettimeofday(&t2, NULL);
		total_time += (float)(t2.tv_sec-t1.tv_sec) + (float)(t2.tv_usec-t1.tv_usec)/1e6;
	}

	void stop() {
		gettimeofday(&t2, NULL);
		total_time += (float)(t2.tv_sec-t1.tv_sec) + (float)(t2.tv_usec-t1.tv_usec)/1e6;
	}

	void clear() {
		total_time = 0;
	}
};
