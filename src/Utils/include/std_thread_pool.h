// thread_pool.h

#pragma once

#include <vector>
#include <queue>
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <future>
#include <memory>
#include <type_traits>
#include <atomic>
#include <sstream>
#include <iostream>

struct TaskArgsGeneral {
    int index_lower{0};
    int index_upper{0};
    int index_job{0};
    void* ptr_any{nullptr};
};

static inline int64_t getTickCountMs()
{
    struct timespec ts;
    clock_gettime ( CLOCK_REALTIME, &ts );
    return static_cast<int64_t>( ( ts.tv_sec * 1e3 ) + ( ts.tv_nsec / 1.e6 ) );
}

class StdThreadPool
{
public:
    explicit StdThreadPool(size_t numThreads = std::thread::hardware_concurrency());
    ~StdThreadPool();
    void begin(const std::string& _name_job);

    template<class F, class... Args>
    // auto push_task(F&& f, Args&&... args) 
    //     -> std::future<typename std::invoke_result<F, Args...>::type>
    void push_task(F&& f, Args&&... args) 
    {
        #if 0
        if(m_name_job.empty() && 0 == m_tasks_queued)
        {
            std::stringstream _stream;
            _stream << "job " << ++m_num_job;
            begin(_stream.str());
        }
        #endif
        
        if(++m_tasks_queued > static_cast<size_t>(m_num_threads))
        {
            throw std::runtime_error("task size must not exceed processor threads...");
        }
        using return_type = typename std::invoke_result<F, Args...>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
            
        // std::future<return_type> res = task->get_future();
        {
            std::lock_guard<std::mutex> lock(m_mutex_queue);
            if(m_stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");

            m_queue_tasks.emplace([this, task](){
                (*task)();
                ++m_tasks_completed;
                m_condition_done.notify_all();
            });
        }
        m_condition_run.notify_one();
        // return res;
    }

    void wait_for_tasks();
    inline const int& get_thread_count()
    {
        return m_num_threads;
    }

    void parallelize(const int& _num_works, const std::function<void(const TaskArgsGeneral& _args)>& _method, void* _ptr_any = nullptr);

protected:
    void thread_worker();
    // void clear_job_done();
    // bool are_all_jobs_done();

private:
    std::vector<std::thread> m_list_thread_workers;
    std::queue<std::function<void()>> m_queue_tasks;
    
    std::mutex m_mutex_queue;
    std::condition_variable m_condition_run;
    std::condition_variable m_condition_done;
    bool m_stop{false};
    std::atomic<size_t> m_tasks_queued{0};
    std::atomic<size_t> m_tasks_completed{0};
    int m_num_threads{0};
    std::string m_name_job;
    // std::vector<uint8_t> m_list_job_done;
    // std::mutex m_mutex_job_done;
    
    static int m_num_job;
};

