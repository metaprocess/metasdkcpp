// std_thread_pool.cpp

#include "std_thread_pool.h"
#include <sys/prctl.h>
#include <sstream>
#include <iostream>
#include "Utils.h"

// int StdThreadPool::m_num_job{0};
StdThreadPool::StdThreadPool(size_t numThreads) : m_stop(false), m_num_threads(numThreads) {
    // m_list_job_done.resize(numThreads);
    // clear_job_done();
    for(size_t i = 0; i < numThreads; ++i)
    {
        m_list_thread_workers.emplace_back(
            [this, i] {
                std::stringstream _ss;
                _ss << "thrdpool." << i ;
                prctl(PR_SET_NAME, _ss.str().c_str(), 0, 0, 0);
                thread_worker();
                // std::cerr << "exiting '" << _ss.str() << "'\n";
            }
        );
    }
}

/**
 * @brief Destructor for StdThreadPool.
 */
StdThreadPool::~StdThreadPool() {
    {
        std::lock_guard<std::mutex> lock(m_mutex_queue);
        m_stop = true;
    }
    m_condition_run.notify_all();
    for (std::thread& worker : m_list_thread_workers) {
        worker.join();
    }
}

/**
 * @brief Initiates a new job with the specified name.
 */
void StdThreadPool::begin(const std::string& _name_job) {
    m_name_job = _name_job;
}

// template<class F, class... Args>
// auto thread_pool::push_task(F&& f, Args&&... args) 
//     -> std::future<typename std::invoke_result<F, Args...>::type>
// {
//     using return_type = typename std::invoke_result<F, Args...>::type;

//     auto task = std::make_shared<std::packaged_task<return_type()>>(
//         std::bind(std::forward<F>(f), std::forward<Args>(args)...)
//     );
        
//     std::future<return_type> res = task->get_future();
//     {
//         std::lock_guard<std::mutex> lock(m_mutex_queue);
//         if(m_stop)
//             throw std::runtime_error("enqueue on stopped ThreadPool");

//         ++m_tasks_total;
//         m_queue_tasks.emplace([this, task](){
//             (*task)();
//             ++m_tasks_completed;
//             if(m_tasks_completed == m_tasks_total) {
//                 m_condition_done.notify_all();
//             }
//         });
//     }
//     m_condition_run.notify_one();
//     return res;
// }

void StdThreadPool::parallelize(const int &_num_works, const std::function<void(const TaskArgsGeneral &_args)> &_method, void *_ptr_any)
{
    // const int _num_works = _num_works;
    m_num_works = std::min(_num_works, m_num_threads);
    int chunk_size = _num_works / m_num_threads;
    if(0 == chunk_size)
    {
        for(int _iter_job = 0; _iter_job < _num_works; _iter_job++)
        {
            TaskArgsGeneral _args{
                .index_lower = _iter_job,
                .index_upper = _iter_job + 1,
                .index_job = _iter_job,
                .ptr_any = _ptr_any
            };
            push_task(_method, _args);
        }
    }
    else
    {
        for(int _iter_job = 0; _iter_job < m_num_threads; _iter_job++)
        {
            TaskArgsGeneral _args{
                .index_lower = chunk_size * _iter_job,
                .index_upper = (chunk_size * (_iter_job+1)) + (_iter_job == (m_num_threads - 1)) * (_num_works - (m_num_threads * chunk_size)),
                .index_job = _iter_job,
                .ptr_any = _ptr_any
            };
            push_task(_method, _args);
        }
    }

    
    // std::cerr << "all tasks pushed\n";
    wait_for_tasks();
}

void StdThreadPool::thread_worker()
{
    while(!m_stop) {
        std::function<void()> task;
        {
            std::unique_lock<std::mutex> lock(m_mutex_queue);
            m_condition_run.wait(
                lock,
                [this]{
                    return m_stop || !m_queue_tasks.empty();
                }
            );
            if(m_queue_tasks.empty())
            {
                continue;
            }
            if(m_stop)
            {
                break;
            }
            task = std::move(m_queue_tasks.front());
            m_queue_tasks.pop();
        }
        task();
        if(m_tasks_completed < static_cast<size_t>(m_num_threads))
        {
            m_condition_run.notify_one();
        }
    }
}

// void thread_pool::clear_job_done()
// {
//     bzero(m_list_job_done.data(), m_list_job_done.size());
// }

// bool thread_pool::are_all_jobs_done()
// {
//     for(const auto& _flag: m_list_job_done)
//     {
//         if(!_flag)
//         {
//             return false;
//         }
//     }
//     return true;
// }

void StdThreadPool::wait_for_tasks() {
    if(m_tasks_completed < static_cast<size_t>(m_num_works))
    {
        std::unique_lock<std::mutex> lock(m_mutex_queue);
        m_condition_done.wait(lock, [this]{
            bool _done = (m_tasks_completed >= static_cast<size_t>(m_num_works));
            // std::cerr << Utils::getTickCountMs() << "::" <<this->m_name_job <<" :: m_tasks_completed = " << m_tasks_completed << "\n";
            return _done;
        });
        // m_condition_done.wait(lock);
    }
    // while(!are_all_jobs_done())
    // {
    //     Utils::microSleep(10);
    // }
    m_tasks_completed = 0;
    m_tasks_queued = 0;
    m_num_works = 0;
    // clear_job_done();
    m_name_job.clear();
}

// inline int StdThreadPool::get_thread_count()
// {
//     return m_num_threads;
// }
