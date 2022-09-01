#ifndef MY_THREADS_H
#define MY_THREADS_H

#include <vector>
#include <numeric>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <chrono>

class work
{
public:
    virtual ~work() {}
    virtual void execute() = 0;
};

typedef work* work_ptr;

class thread_pool
{
protected:
    
    bool single_threaded;
    
    std::atomic<int> cur, done;
    std::vector<work_ptr> todo;
    
    std::mutex work_mutex;
    std::condition_variable work_available;
    
    std::mutex done_mutex;
    std::condition_variable work_done;
    
    bool destroy;
    
    class thread_worker 
    {
    protected:
        thread_pool *pool;
        int id;
        std::thread thread;
        
    public:
        thread_worker(thread_pool *pool, int id) : pool(pool), id(id) {
            thread = std::thread(&thread_pool::thread_worker::execute, this);
        };
        ~thread_worker() { };
        
        void execute() {
            
            while (true) {
                
                {
                    std::unique_lock<std::mutex> lk(pool->work_mutex);
                    pool->work_available.wait(lk, [&]{ return pool->cur < (int)pool->todo.size() || pool->destroy; });
                }
                
                if (pool->destroy)
                    return;
                
                int work_cnt = pool->todo.size();
                while (true) {
                    int next = pool->cur++;
                    if (next < work_cnt) {
                        work_ptr w = pool->todo[next];
                        w->execute();
                    } else {
                        break;
                    }
                    
                    pool->done_mutex.lock();
                    int done = ++(pool->done);
                    pool->done_mutex.unlock();
                    
                    if (done == work_cnt) {
                        pool->work_done.notify_one();
                    }
                }
            }
        }
        
        void join() {
            thread.join();
        }
    };
    
    std::vector<thread_worker*> workers;
    
public:
    
    thread_pool(int num_threads) : single_threaded(false), cur(0), done(0), todo(0), destroy(false) {
        if (num_threads < 1)
            throw new std::invalid_argument("number of threads has to be at least 1");
        
        if (num_threads == 1) {
            single_threaded = true;
        } else {
            for (int i = 0; i < num_threads; i++)
                workers.push_back(new thread_worker(this, i));
        }
    }
    
    ~thread_pool() {
        destroy = true;
        work_available.notify_all();
        
        for (auto w : workers) {
            w->join();
            delete w;
        }
    }
    
    
    void do_work(const std::vector<work_ptr> &todo) {
        
        if (single_threaded) {
            for (auto w : todo)
                w->execute();
        } else {
            
            {
                std::unique_lock<std::mutex> lk(work_mutex);
                this->todo = todo;
                cur = 0;
                done = 0;
            }
            
            work_available.notify_all();
            
            {
                std::unique_lock<std::mutex> lk(done_mutex);
                work_done.wait(lk, [&]{ return done == (int)todo.size(); });
            }
            
        }
    }
};


#endif

