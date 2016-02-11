/*
 * ThreadPool.cpp
 *
 * Copyright (c) 2016 Tobias Wood
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Adapted from http://alexagafonov.com/2015/05/05/thread-pool-implementation-in-c-11/
 */
 
 #include "ThreadPool.h"
 
 namespace QI {
     
    ThreadPool::ThreadPool(const size_t nThreads, const bool d) : m_stopping(false), m_debug(d) {
        for (size_t i = 0; i < nThreads; i++) {
            m_threads.emplace_back(std::thread(&ThreadPool::invokeThread, this));
            if (m_debug) std::cout << "Emplaced thread " << m_threads.back().get_id() << std::endl;
        }
        if (m_debug) std::cout << "Constructed thread pool with " << nThreads << " threads" << std::endl;
    }
    
    ThreadPool::~ThreadPool() {
        if (m_debug) std::cout << "Destructing thread pool" << std::endl;
        { // Lock will exist in this scope
            std::unique_lock<std::mutex> lock(m_tasksMutex);
            m_stopping = true;
        }
        m_threadCondition.notify_all();
        for (std::thread &t : m_threads) {
            if (m_debug) std::cout << "Joining thread " << t.get_id() << std::endl;
            t.join();
        }
        m_threads.empty();
        if (m_debug) std::cout << "Destructed thread pool" << std::endl;
    }
    
    void ThreadPool::setDebug(const bool d) { m_debug = d; }
    
    void ThreadPool::enqueue(TFunc f) {
        {   // Wait for space in the queue
            std::unique_lock<std::mutex> qLock(m_queueMutex);
            m_queueCondition.wait(qLock, [this]{ return m_tasks.size() < 2 * m_threads.size(); });
        }
        
        { // Lock will exist in this scope
            std::unique_lock<std::mutex> tLock(m_tasksMutex);
            m_tasks.push(f);
            if (m_debug) std::cout << "Pushed task" << &f << std::endl;
        }
        m_threadCondition.notify_one(); // Wake up a thread
    }
    
    void ThreadPool::invokeThread() {
        TFunc task;
        if (m_debug) std::cout << "Thread started" << std::endl;
        while (true) {
            { // Lock will exist in this scope
                std::unique_lock<std::mutex> lock(m_tasksMutex);
                m_threadCondition.wait(lock, [this]{ return !m_tasks.empty() || m_stopping; });
                if (m_stopping && m_tasks.empty()) {
                    return;
                }
                task = m_tasks.front();
                m_tasks.pop();
                m_queueCondition.notify_one(); // Wake up main thread if something is waiting to enqueue
                if (m_debug) std::cout << "Starting task " << &task << std::endl;
            }
            task();
            
        }
        if (m_debug) std::cout << "Thread stopped" << std::endl;
    }
    
 } // End namespace QI