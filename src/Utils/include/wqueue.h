/*
   wqueue.h

   Worker thread queue based on the Standard C++ library list
   template class.
*/

#ifndef __wqueue_h__
#define __wqueue_h__

#include <pthread.h>
#include <list>
#include <unistd.h>
#include <iostream>
#ifndef nullptr
#define nullptr NULL
#endif

//using namespace std;

template <typename T> class wqueue
{

  public:
    wqueue() {
        pthread_mutex_init(&m_mutex, nullptr);
        pthread_cond_init(&m_condv, nullptr);
		m_bRunning = true;
    }
    ~wqueue() {
		m_bRunning = false;
		pthread_cond_signal(&m_condv);
		usleep(5000);
        pthread_mutex_destroy(&m_mutex);
        pthread_cond_destroy(&m_condv);
//		std::cerr << "destructor " << __FUNCTION__ << "\n";
    }
    void add(T item) {
        pthread_mutex_lock(&m_mutex);
        m_queue.push_back(item);
        pthread_cond_signal(&m_condv);
        pthread_mutex_unlock(&m_mutex);
    }
    T remove(const bool* parentRunning = nullptr) {
		T item = nullptr;
        pthread_mutex_lock(&m_mutex);
		while(true)
		{
			while (m_queue.size() == 0 ) 
			{
				pthread_cond_wait(&m_condv, &m_mutex);
				if(!m_bRunning)
				{
					std::cerr << "return from wqueue " << __FUNCTION__ << "\n";
					break;
				}
			}
			if(m_queue.size() > 0)
			{
				item = m_queue.front();
				if(!item )
				{
					if( !parentRunning || !(*parentRunning) )
						m_queue.pop_front();
					else 
					{
						usleep(1000);
						pthread_cond_signal(&m_condv);
						pthread_mutex_unlock(&m_mutex);
						continue;
					}
				}
				else
				{
					m_queue.pop_front();
				}
				
			}
			else
			{
				continue;
				pthread_mutex_unlock(&m_mutex);
			}
			
			break;
		}
		pthread_mutex_unlock(&m_mutex);
		return item;
    }
    T pop_first()
    {

        T _t= nullptr;
        pthread_mutex_lock(&m_mutex);
        _t = m_queue.front();
        m_queue.pop_front();
        pthread_mutex_unlock(&m_mutex);
        return _t;
    }
    int size() {
        pthread_mutex_lock(&m_mutex);
        int size = m_queue.size();
        pthread_mutex_unlock(&m_mutex);
        return size;
    }
    
    void stop() { m_bRunning = false;}
    
private:
	std::list<T>          m_queue;
    pthread_mutex_t  m_mutex;
    pthread_cond_t   m_condv; 
	bool m_bRunning;
};

template <typename T> class wqueue2
{

  public:
	wqueue2() {
		pthread_mutex_init(&m_mutex, nullptr);
		pthread_cond_init(&m_condv, nullptr);
		m_bRunning = true;
	}
	~wqueue2() {
		m_bRunning = false;
		pthread_cond_signal(&m_condv);
		usleep(5000);
		pthread_mutex_destroy(&m_mutex);
		pthread_cond_destroy(&m_condv);
//		std::cerr << "destructor " << __FUNCTION__ << "\n";
	}
	void add(T item) {
		pthread_mutex_lock(&m_mutex);
		m_queue.push_back(item);
		pthread_cond_signal(&m_condv);
		pthread_mutex_unlock(&m_mutex);
	}
	T remove(const bool* parentRunning = nullptr) {
		T item;// = nullptr;
		pthread_mutex_lock(&m_mutex);
		while(true)
		{
			while (m_queue.size() == 0 )
			{
				pthread_cond_wait(&m_condv, &m_mutex);
				if(!m_bRunning)
				{
					std::cerr << "return from wqueue " << __FUNCTION__ << "\n";
					break;
				}
			}
			if(m_queue.size() > 0)
			{
				item = m_queue.front();
//				if(!item )
				{
					if( !parentRunning || !(*parentRunning) )
						m_queue.pop_front();
					else
					{
						usleep(1000);
						pthread_cond_signal(&m_condv);
						pthread_mutex_unlock(&m_mutex);
						continue;
					}
				}
//				else
//				{
//					m_queue.pop_front();
//				}

			}
			else
			{
				continue;
				pthread_mutex_unlock(&m_mutex);
			}

			break;
		}
		pthread_mutex_unlock(&m_mutex);
		return item;
	}
    T pop_first()
    {

        pthread_mutex_lock(&m_mutex);
        T _t = m_queue.front();
        m_queue.pop_front();
        pthread_mutex_unlock(&m_mutex);
        return _t;
    }
	int size() {
		pthread_mutex_lock(&m_mutex);
		int size = m_queue.size();
		pthread_mutex_unlock(&m_mutex);
		return size;
	}

	void stop() { m_bRunning = false;}

private:
	std::list<T>          m_queue;
	pthread_mutex_t  m_mutex;
	pthread_cond_t   m_condv;
	bool m_bRunning;
};

#endif
