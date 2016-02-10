/*
 * ThreadPool.h
 *
 * Copyright (c) 2016 Tobias Wood
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Adapted from http://alexagafonov.com/2015/05/05/thread-pool-implementation-in-c-11/
 */
 
#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <vector>
#include <queue>

namespace QI {
    
class ThreadPool {
    
private:
    std::vector<std::thread> m_threads;
    std::queue<std::function<void ()>> m_tasks;
    std::mutex m_tasksMutex, m_queueMutex;
    std::condition_variable m_queueCondition, m_threadCondition;
    
    bool m_stopping, m_debug;
    
    void invokeThread();
    
public:
    typedef std::function<void (void)> TFunc;
    
    ThreadPool(const size_t nThreads);
    ~ThreadPool();
    
    void enqueue(TFunc f);
    void setDebug(const bool d);
};
    
    
} // End namespace QI

#endif // End THREAD_POOL_H