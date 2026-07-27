#ifndef QONMODULE_TIMER_H
#define QONMODULE_TIMER_H

#include <cstdint>

class HighResTimer
{
 public:
  HighResTimer() = default;
  HighResTimer(bool start);
  ~HighResTimer() = default;
  void Start();
  void Stop();
  void Abort();
  void Reset();
  void ResetStart();
  double GetElapsedTime();
  double GetCurrentElapsedTime(bool reset = false);
  void StopAndStart(HighResTimer& startTimer);
  int32_t IsRunning() { return running; }
  void AddTime(double t);

 private:
  double ElapsedTime = 0.;
  double StartTime = 0.;
  int32_t running = 0;

  static double GetFrequency();
  static double GetTime();
  static double Frequency;
};

#endif
