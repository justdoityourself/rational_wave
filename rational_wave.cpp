/* Copyright (C) 2020 D8DATAWORKS - All Rights Reserved */

#ifdef TEST_RUNNER


#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "rw/test.hpp"

int main(int argc, char* argv[])
{
    return Catch::Session().run(argc, argv);
}


#endif //TEST_RUNNER

#ifdef BENCHMARK_RUNNER


#define PICOBENCH_IMPLEMENT
#include "picobench.hpp"
#include "rw/benchmark.hpp"

int main(int argc, char* argv[])
{
    picobench::runner runner;
    return runner.run();
}


#endif //BENCHMARK_RUNNER



#if ! defined(BENCHMARK_RUNNER) && ! defined(TEST_RUNNER) && ! defined(OTHER)

int main(int argc, char* argv[])
{
    return 0;
}


#endif //! defined(BENCHMARK_RUNNER) && ! defined(TEST_RUNNER)