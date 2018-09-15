#ifndef PARALLEL_FOR_H
#define PARALLEL_FOR_H
#if defined (OPENMP_THREADING)
#define PARALLEL_FOR(nthreads,LOOP_END,O) {          		\
	_Pragma("omp parallel num_threads(nthreads)")			\
	{							\
	_Pragma("omp for")						\
	for (int loop_i=0; loop_i<LOOP_END; loop_i++) {		\
		O;						\
	}							\
	}							\
}

#elif defined (NATIVE_THREADING)
#define PARALLEL_FOR(nthreads,LOOP_END,O) {          \
        std::vector<std::thread> threads(nthreads);			            \
        for (int t = 0; t < nthreads; t++) {					\
            threads[t] = std::thread(std::bind(					\
                    [&](const int bi, const int ei, const int t) { 		\
		for(int loop_i = bi;loop_i<ei;loop_i++) {   O;  }	\
				    },t*LOOP_END/nthreads,(t+1)==nthreads?LOOP_END:(t+1)*LOOP_END/nthreads,t)); 						\
	}								\
        std::for_each(threads.begin(),threads.end(),[](std::thread& x){x.join();});\
}

#else 
#define PARALLEL_FOR(nthreads,LOOP_END,O) {          \
	for (int loop_i=0; loop_i<LOOP_END; loop_i++) {		\
		O;				\
	}			\
}
#endif
#endif
