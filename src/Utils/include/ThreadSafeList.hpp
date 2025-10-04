#include <list>
#include <mutex>
#include <shared_mutex>
#include <condition_variable>

template <typename T>
class ThreadSafeList {
private:
    std::list<T> list_;
    mutable std::shared_mutex mutex_;
    std::condition_variable_any cv_;

public:
    std::vector<T> get_all_elements() const {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        return std::vector<T>(list_.begin(), list_.end());
    }
    // Thread-safe iterator class
    class Iterator {
    private:
        ThreadSafeList<T>& list_;
        typename std::list<T>::iterator it_;
        std::unique_ptr<std::shared_lock<std::shared_mutex>> lock_;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T*;
        using reference = T&;

        // Constructor
        Iterator(ThreadSafeList<T>& list, typename std::list<T>::iterator it)
            : list_(list), it_(it), 
              lock_(std::make_unique<std::shared_lock<std::shared_mutex>>(list.mutex_)) {}

        // Move constructor
        Iterator(Iterator&& other) noexcept
            : list_(other.list_),
              it_(std::move(other.it_)),
              lock_(std::move(other.lock_)) {}

        // Move assignment operator
        Iterator& operator=(Iterator&& other) noexcept {
            if (this != &other) {
                list_ = other.list_;
                it_ = std::move(other.it_);
                lock_ = std::move(other.lock_);
            }
            return *this;
        }

        // Delete copy constructor and assignment
        Iterator(const Iterator&) = delete;
        Iterator& operator=(const Iterator&) = delete;

        Iterator& operator++() {
            ++it_;
            return *this;
        }

        Iterator operator++(int) {
            Iterator tmp(std::move(*this));
            ++it_;
            return tmp;
        }

        reference operator*() { return *it_; }
        pointer operator->() { return &(*it_); }

        bool operator==(const Iterator& other) const {
            return it_ == other.it_;
        }

        bool operator!=(const Iterator& other) const {
            return !(*this == other);
        }
    };

    ThreadSafeList() = default;

    // Begin iterator
    Iterator begin() {
        return Iterator(*this, list_.begin());
    }

    // End iterator
    Iterator end() {
        return Iterator(*this, list_.end());
    }

    // Rest of the methods remain the same...
    void push_back(const T& element) {
        std::unique_lock<std::shared_mutex> lock(mutex_);
        list_.push_back(element);
        cv_.notify_all();
    }

    void push_front(const T& element) {
        std::unique_lock<std::shared_mutex> lock(mutex_);
        list_.push_front(element);
        cv_.notify_all();
    }

    bool pop_front(T& element) {
        std::unique_lock<std::shared_mutex> lock(mutex_);
        if (list_.empty()) {
            return false;
        }
        element = list_.front();
        list_.pop_front();
        return true;
    }

    bool pop_back(T& element) {
        std::unique_lock<std::shared_mutex> lock(mutex_);
        if (list_.empty()) {
            return false;
        }
        element = list_.back();
        list_.pop_back();
        return true;
    }

    bool peek_front(T& element) const {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        if (list_.empty()) {
            return false;
        }
        element = list_.front();
        return true;
    }

    bool peek_back(T& element) const {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        if (list_.empty()) {
            return false;
        }
        element = list_.back();
        return true;
    }

    bool empty() const {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        return list_.empty();
    }

    size_t size() const {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        return list_.size();
    }

    void wait_until_not_empty() const {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        cv_.wait(lock, [this] { return !list_.empty(); });
    }
};