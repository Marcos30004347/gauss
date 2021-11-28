#ifndef TEST_H
#define TEST_H

#include <assert.h>
#include <chrono>
#include <iostream>

std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;
#define TIMER_START(test_name) begin = std::chrono::steady_clock::now();
#define TIMER_STOP(test_name)                                                  \
  end = std::chrono::steady_clock::now();                                      \
  std::cout << test_name << " PASS in "                                        \
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -     \
                                                                     begin)    \
                   .count()                                                    \
            << "ms" << std::endl;

#define TEST(test_name)                                                        \
  TIMER_START(#test_name) test_name();                                         \
  TIMER_STOP(#test_name)

#endif
