/*
 * Copyright (c) 2010-2016 Richard Geldreich, Jr. and Binomial LLC
 * Copyright (c) 2020 FrozenStorm Interactive, Yoann Potinet
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any damages
 * arising from the use of this software.
 * 
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 * 
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation or credits
 *    is required.
 * 
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 
 * 3. This notice may not be removed or altered from any source distribution.
 */

#pragma once

#include "crn_core.h"

#if CRNLIB_USE_PTHREADS_API

#include "crn_atomics.h"

#if CRNLIB_NO_ATOMICS
#error No atomic operations defined in crn_platform.h!
#endif

#if defined(CRN_OS_DARWIN)
#include <os/lock.h>
#endif

#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>

#include "crn_export.h"

namespace crnlib
{
    // g_number_of_processors defaults to 1. Will be higher on multicore machines.
    CRN_EXPORT extern uint g_number_of_processors;

    CRN_EXPORT void crn_threading_init();

    typedef uint64 crn_thread_id_t;
    CRN_EXPORT crn_thread_id_t crn_get_current_thread_id();

    CRN_EXPORT void crn_sleep(unsigned int milliseconds);

    CRN_EXPORT uint crn_get_max_helper_threads();

    class CRN_EXPORT mutex
    {
        mutex(const mutex&);
        mutex& operator=(const mutex&);
    public:
        mutex(unsigned int spin_count = 0);
        ~mutex();
        void lock();
        void unlock();
        void set_spin_count(unsigned int count);

    private:
        pthread_mutex_t m_mutex;

#ifdef CRNLIB_BUILD_DEBUG
        unsigned int m_lock_count;
#endif
    };

    class CRN_EXPORT scoped_mutex
    {
        scoped_mutex(const scoped_mutex&);
        scoped_mutex& operator=(const scoped_mutex&);
    public:
        inline scoped_mutex(mutex& m):
            m_mutex(m)
        {
            m_mutex.lock();
        }
        inline ~scoped_mutex()
        {
            m_mutex.unlock();
        }

    private:
        mutex& m_mutex;
    };

    class CRN_EXPORT semaphore
    {
        CRNLIB_NO_COPY_OR_ASSIGNMENT_OP(semaphore);
    public:
        semaphore(long initialCount = 0, long maximumCount = 1, const char* pName = nullptr);
        ~semaphore();

        void release(long releaseCount = 1);
        void try_release(long releaseCount = 1);
        bool wait(uint32 milliseconds = cUINT32_MAX);

    private:
        sem_t m_sem;
    };

    class CRN_EXPORT spinlock
    {
    public:
        spinlock();
        ~spinlock();

        void lock();
        void unlock();

    private:
#if defined(CRN_OS_LINUX)
        pthread_spinlock_t m_spinlock;
#elif defined(CRN_OS_DARWIN)
        os_unfair_lock m_spinlock;
#else
        int m_spinlock;
#endif
    };

    class CRN_EXPORT scoped_spinlock
    {
        scoped_spinlock(const scoped_spinlock&);
        scoped_spinlock& operator=(const scoped_spinlock&);
    public:
        inline scoped_spinlock(spinlock& lock):
            m_lock(lock)
        {
            m_lock.lock();
        }
        inline ~scoped_spinlock()
        {
            m_lock.unlock();
        }

    private:
        spinlock& m_lock;
    };

    template <typename T, uint cMaxSize>
    class tsstack
    {
    public:
        inline tsstack():
            m_top(0)
        {
        }

        inline ~tsstack()
        {
        }

        inline void clear()
        {
            m_spinlock.lock();
            m_top = 0;
            m_spinlock.unlock();
        }

        inline bool try_push(const T& obj)
        {
            bool result = false;
            m_spinlock.lock();
            if (m_top < (int)cMaxSize)
            {
                m_stack[m_top++] = obj;
                result = true;
            }
            m_spinlock.unlock();
            return result;
        }

        inline bool pop(T& obj)
        {
            bool result = false;
            m_spinlock.lock();
            if (m_top > 0)
            {
                obj = m_stack[--m_top];
                result = true;
            }
            m_spinlock.unlock();
            return result;
        }

    private:
        spinlock m_spinlock;
        T m_stack[cMaxSize];
        int m_top;
    };

    class CRN_EXPORT task_pool
    {
    public:
        task_pool();
        task_pool(uint num_threads);
        ~task_pool();

        enum { cMaxThreads = 16 };
        bool init(uint num_threads);
        void deinit();

        inline uint get_num_threads() const
        {
            return m_num_threads;
        }
        inline uint32 get_num_outstanding_tasks() const
        {
            return m_total_submitted_tasks - m_total_completed_tasks;
        }

        // C-style task callback
        typedef void (*task_callback_func)(uint64 data, void* pData_ptr);
        bool queue_task(task_callback_func pFunc, uint64 data = 0, void* pData_ptr = nullptr);

        class executable_task
        {
        public:
            virtual void execute_task(uint64 data, void* pData_ptr) = 0;
        };

        // It's the caller's responsibility to delete pObj within the execute_task() method, if needed!
        bool queue_task(executable_task* pObj, uint64 data = 0, void* pData_ptr = nullptr);

        template <typename S, typename T>
        inline bool queue_object_task(S* pObject, T pObject_method, uint64 data = 0, void* pData_ptr = nullptr);

        template <typename S, typename T>
        inline bool queue_multiple_object_tasks(S* pObject, T pObject_method, uint64 first_data, uint num_tasks, void* pData_ptr = nullptr);

        void join();

    private:
        struct task
        {
            inline task():
                m_data(0),
                m_pData_ptr(nullptr),
                m_pObj(nullptr),
                m_flags(0)
            {
            }

            uint64 m_data;
            void* m_pData_ptr;

            union {
                task_callback_func m_callback;
                executable_task* m_pObj;
            };

            uint m_flags;
        };

        tsstack<task, cMaxThreads> m_task_stack;

        uint m_num_threads;
        pthread_t m_threads[cMaxThreads];

        // Signalled whenever a task is queued up.
        semaphore m_tasks_available;

        // Signalled when all outstanding tasks are completed.
        semaphore m_all_tasks_completed;

        enum task_flags
        {
            cTaskFlagObject = 1
        };

        volatile atomic32_t m_total_submitted_tasks;
        volatile atomic32_t m_total_completed_tasks;
        volatile atomic32_t m_exit_flag;

        void process_task(task& tsk);

        static void* thread_func(void* pContext);
    };

    enum object_task_flags
    {
        cObjectTaskFlagDefault = 0,
        cObjectTaskFlagDeleteAfterExecution = 1
    };

    template <typename T>
    class object_task : public task_pool::executable_task
    {
    public:
        object_task(uint flags = cObjectTaskFlagDefault):
            m_pObject(nullptr),
            m_pMethod(nullptr),
            m_flags(flags)
        {
        }

        typedef void (T::* object_method_ptr)(uint64 data, void* pData_ptr);

        object_task(T* pObject, object_method_ptr pMethod, uint flags = cObjectTaskFlagDefault):
            m_pObject(pObject),
            m_pMethod(pMethod),
            m_flags(flags)
        {
            CRNLIB_ASSERT(pObject && pMethod);
        }

        void init(T* pObject, object_method_ptr pMethod, uint flags = cObjectTaskFlagDefault)
        {
            CRNLIB_ASSERT(pObject && pMethod);

            m_pObject = pObject;
            m_pMethod = pMethod;
            m_flags = flags;
        }

        T* get_object() const
        {
            return m_pObject;
        }
        object_method_ptr get_method() const
        {
            return m_pMethod;
        }

        virtual void execute_task(uint64 data, void* pData_ptr)
        {
            (m_pObject->*m_pMethod)(data, pData_ptr);

            if (m_flags & cObjectTaskFlagDeleteAfterExecution)
            {
                crnlib_delete(this);
            }
        }

    protected:
        T* m_pObject;

        object_method_ptr m_pMethod;

        uint m_flags;
    };

    template <typename S, typename T>
    inline bool task_pool::queue_object_task(S* pObject, T pObject_method, uint64 data, void* pData_ptr)
    {
        object_task<S>* pTask = crnlib_new<object_task<S> >(pObject, pObject_method, cObjectTaskFlagDeleteAfterExecution);
        if (!pTask)
        {
            return false;
        }
        return queue_task(pTask, data, pData_ptr);
    }

    template <typename S, typename T>
    inline bool task_pool::queue_multiple_object_tasks(S* pObject, T pObject_method, uint64 first_data, uint num_tasks, void* pData_ptr)
    {
        CRNLIB_ASSERT(pObject);
        CRNLIB_ASSERT(num_tasks);
        if (!num_tasks)
        {
            return true;
        }

        bool status = true;

        uint i;
        for (i = 0; i < num_tasks; i++)
        {
            task tsk;

            tsk.m_pObj = crnlib_new<object_task<S> >(pObject, pObject_method, cObjectTaskFlagDeleteAfterExecution);
            if (!tsk.m_pObj)
            {
                status = false;
                break;
            }

            tsk.m_data = first_data + i;
            tsk.m_pData_ptr = pData_ptr;
            tsk.m_flags = cTaskFlagObject;

            atomic_increment32(&m_total_submitted_tasks);

            if (!m_task_stack.try_push(tsk))
            {
                atomic_increment32(&m_total_completed_tasks);

                status = false;
                break;
            }
        }

        if (i)
        {
            m_tasks_available.release(i);
        }

        return status;
    }
}  // namespace crnlib

#endif  // CRNLIB_USE_PTHREADS_API
