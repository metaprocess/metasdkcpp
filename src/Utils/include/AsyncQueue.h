#ifndef QASYNCQUEUE_H
#define QASYNCQUEUE_H


#include <QMutex>
#include <QQueue>
#include <QtCore>
#include "Definitions.h"
#ifdef __linux
	#include <unistd.h>
#endif
template <class T>
	/**
	 * @brief this is a container class that containes member object of type Class T . it's algorithm is FIFO "first in first out" 
	 * 
	 */
class AsyncQueue
{
public:

	AsyncQueue(int32_t maxSize = 10000)
	{
		m_MaxSize = maxSize;
		
		m_Queue = new QQueue<T>();
		m_Mutex = new QMutex();
	}
	
	void SetMaxSize(int maxSize)
	{
		m_MaxSize = maxSize;
	}
	
	~AsyncQueue()
	{
		clear();
		delete m_Mutex;
		delete m_Queue;
	}
	
	void TakeFirst(T& data)
	{
		QMutexLocker lck(m_Mutex);
		// 		if( m_Queue->count() > 0 )
		// 		m_Queue->
		data = m_Queue->takeFirst();
		
		// 		return b;
	}
	
	void TakeLast(T& data)
	{
		QMutexLocker lck(m_Mutex);
		// 		if( m_Queue->count() > 0 )
		// 		m_Queue->
		data = m_Queue->takeLast();
		
		// 		return b;
	}
	void Push(const T& data)
	{
#ifdef __linux
		while(true)
		{
			if(IsFull())
				sleep(1);
			else
				break;
		}
#endif
		m_Mutex->lock();
		m_Queue->enqueue(data);
		m_Mutex->unlock();
	}
	
	int32_t maxSize(void)
	{
		return m_MaxSize;
	}
	
	void clear()
	{
		m_Mutex->lock();
		m_Queue->clear();
		m_Mutex->unlock();
	}
	
	void Pull(T& data)
	{
		#ifdef __linux
		while(true)
		{
			if(IsEmpty())
				sleep(1);
			else
				break;
		}
#endif
		m_Mutex->lock();
		data = m_Queue->dequeue();
		m_Mutex->unlock();
	}
	
	bool IsFull()
	{
		if(m_MaxSize < 1)
			return true;
		
		m_Mutex->lock();
		int count = m_Queue->count();
		m_Mutex->unlock();
        return (count >= m_MaxSize);
	}
	
	bool IsEmpty()
	{
		bool empty = m_Queue->isEmpty();
		return empty;
	}
	
	int32_t Count()
	{
		int32_t cnt;
		m_Mutex->lock();
		cnt = m_Queue->count();
		m_Mutex->unlock();
		return cnt;
	}
	
private:
	QMutex* m_Mutex;
	QQueue<T>* m_Queue;
	int32_t m_MaxSize;
	
};

#endif // QASYNCQUEUE_H
