#!/usr/bin/env python
import multiprocessing, os, signal, time, Queue, functions

def do_work(value):
    pid = os.getpid()
    cmd = " echo \"%d part 1\"; sleep 3; echo \"%d part 2\"; sleep 3; echo \"%d complete\" " % (value, value, value)
    result = functions.srun(cmd)
    print("%d was %s" % (value, str(result)))

def worker_function(job_queue):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    while not job_queue.empty():
        try:
            job = job_queue.get(block=False)
            do_work(job)
        except Queue.Empty:
            pass
        #except KeyboardInterrupt: pass

def main():
    job_queue = multiprocessing.Queue()

    for i in range(40):
        job_queue.put(i)

    workers = []
    for i in range(3):
        tmp = multiprocessing.Process(target=worker_function,
                                      args=(job_queue,))
        tmp.start()
        workers.append(tmp)

    try:
        for worker in workers:
            worker.join()
    except KeyboardInterrupt:
        print 'parent received ctrl-c; continuing jobs already submitted to slurm'
        for worker in workers:
            worker.terminate()
            worker.join()

if __name__ == "__main__":
    main()