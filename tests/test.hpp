#ifndef TEST_H
#define TEST_H

#include <chrono>
#include <iostream>
#include <assert.h>

std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;

#define TIMER_START(test_name)                                                 \
  std::cout << test_name << " PASS in ";                                       \
  begin = std::chrono::steady_clock::now();
#define TIMER_STOP                                                             \
  end = std::chrono::steady_clock::now();                                      \
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end -     \
                                                                     begin)    \
                   .count()                                                    \
            << "ms" << std::endl;

#define TEST(test_name)                                                        \
  TIMER_START(#test_name) test_name();                                         \
  TIMER_STOP

#endif
