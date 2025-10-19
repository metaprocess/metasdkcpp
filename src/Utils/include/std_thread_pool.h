// std_thread_pool.h

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

/**
 * @file std_thread_pool.h
 * @brief Header file for the StdThreadPool class, providing a thread pool implementation for concurrent task execution.
 *
 * This file defines the StdThreadPool class, which manages a pool of worker threads to execute tasks concurrently.
 * It leverages standard C++ concurrency features such as std::thread, std::mutex, and std::condition_variable to provide
 * efficient task scheduling and execution. The thread pool supports task queuing, synchronization, and work parallelization.
 */

/**
 * @struct TaskArgsGeneral
 * @brief A general-purpose structure to hold arguments for tasks executed by the thread pool.
 *
 * This structure is used to pass parameters to tasks, particularly for parallelized workloads, providing index ranges,
 * a job identifier, and an optional pointer to arbitrary data.
 */
struct TaskArgsGeneral {
    int index_lower{0};     ///< The lower index of the range for the task.
    int index_upper{0};     ///< The upper index of the range for the task.
    int index_job{0};       ///< The job index associated with the task.
    void* ptr_any{nullptr}; ///< A void pointer to hold any additional data required by the task.
};

/**
 * @brief Retrieves the current time in milliseconds.
 *
 * This function uses the POSIX clock_gettime function to obtain the current time with millisecond precision.
 *
 * @return The current time in milliseconds as an int64_t.
 */
static inline int64_t getTickCountMs() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return static_cast<int64_t>((ts.tv_sec * 1e3) + (ts.tv_nsec / 1.e6));
}

/**
 * @class StdThreadPool
 * @brief A thread pool class for managing concurrent task execution.
 *
 * The StdThreadPool class provides a mechanism to execute tasks concurrently using a fixed number of worker threads.
 * It allows tasks to be pushed to a queue for execution, waits for task completion, and supports parallelizing work
 * across multiple threads. The class uses synchronization primitives to ensure thread safety and efficient resource use.
 */
class StdThreadPool {
public:
    /**
     * @brief Constructs a StdThreadPool with a specified number of threads.
     *
     * Initializes the thread pool with the specified number of worker threads. If no value is provided,
     * it defaults to the number of hardware threads available on the system (std::thread::hardware_concurrency()).
     *
     * @param numThreads The number of worker threads to create in the pool.
     */
    explicit StdThreadPool(size_t numThreads = std::thread::hardware_concurrency());

    /**
     * @brief Destructor for StdThreadPool.
     *
     * Terminates all worker threads by setting the stop flag and joining each thread,
     * ensuring proper cleanup of resources before the object is destroyed.
     */
    ~StdThreadPool();

    /**
     * @brief Initiates a new job with the specified name.
     *
     * Starts a new job by setting the job name, which can be used for tracking or debugging purposes.
     * This method may reset internal states related to the current job.
     *
     * @param _name_job The name of the job to begin.
     */
    void begin(const std::string& _name_job);

    /**
     * @brief Pushes a task to the thread pool for execution.
     *
     * This template method enqueues a callable (e.g., function, lambda) along with its arguments into the task queue.
     * The task is executed by an available worker thread. The method enforces a limit where the number of queued tasks
     * must not exceed the number of threads, throwing an exception if violated.
     *
     * @tparam F The type of the callable.
     * @tparam Args The types of the arguments to the callable.
     * @param f The callable to be executed.
     * @param args The arguments to pass to the callable.
     * @throws std::runtime_error If the task queue exceeds the thread count or if enqueuing on a stopped pool.
     */
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
    /**
     * @brief Retrieves the number of threads in the pool.
     *
     * Provides access to the number of worker threads currently managed by the pool.
     *
     * @return A constant reference to the number of threads.
     */
    inline const int& get_thread_count()
    {
        return m_num_threads;
    }

    /**
     * @brief Distributes a task across multiple threads for parallel execution.
     *
     * Parallelizes a workload by dividing it into a specified number of work items, each processed by a task
     * that invokes the provided method with TaskArgsGeneral arguments. An optional pointer to additional data
     * can be passed to the method.
     *
     * @param _num_works The total number of work items to process.
     * @param _method The method to execute for each work item, accepting a TaskArgsGeneral argument.
     * @param _ptr_any An optional pointer to additional data for the method (default is nullptr).
     */
    void parallelize(const int& _num_works, const std::function<void(const TaskArgsGeneral& _args)>& _method, void* _ptr_any = nullptr);

protected:
    /**
     * @brief The worker function executed by each thread in the pool.
     *
     * Runs in a loop within each worker thread, waiting for tasks to be available in the queue.
     * When a task is dequeued, it is executed, and the loop continues until the stop flag is set.
     */
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
    int m_num_works{0};
    std::string m_name_job;
    // std::vector<uint8_t> m_list_job_done;
    // std::mutex m_mutex_job_done;
    
    // static int m_num_job;
};

