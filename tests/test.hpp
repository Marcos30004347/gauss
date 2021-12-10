#ifndef TEST_H
#define TEST_H

#include <assert.h>
#include <chrono>
#include <iostream>

#define TIMED_SECTION_START(name)                                              \
  {                                                                            \
    std::cout << "timed section " << name << " started\n";                     \
    std::chrono::steady_clock::time_point begin =                              \
        std::chrono::steady_clock::now();

#define TIMED_SECTION_STOP(name)                                               \
  std::chrono::steady_clock::time_point end =                                  \
      std::chrono::steady_clock::now();                                        \
  std::cout << "timed section " << name << " ended in "                        \
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end -      \
                                                                    begin)     \
                   .count()                                                    \
            << "ns" << std::endl;                                              \
  }

#define TEST_START(test_name)                                                  \
  std::chrono::steady_clock::time_point begin =                                \
      std::chrono::steady_clock::now();
#define TEST_STOP(test_name)                                                   \
  std::chrono::steady_clock::time_point end =                                  \
      std::chrono::steady_clock::now();                                        \
  std::cout << test_name << " PASS in "                                        \
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -     \
                                                                     begin)    \
                   .count()                                                    \
            << "ms" << std::endl;

#define TEST(test_name)                                                        \
  {                                                                            \
    TEST_START(#test_name) test_name();                                        \
    TEST_STOP(#test_name)                                                      \
  }

#endif
