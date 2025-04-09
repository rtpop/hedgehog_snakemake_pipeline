import subprocess
import resource
import sys
import time

def measure_resources(command):
    start_time = time.time()
    usage_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    
    result = subprocess.run(command, shell=True)
    
    usage_end = resource.getrusage(resource.RUSAGE_CHILDREN)
    end_time = time.time()
    
    elapsed_time = end_time - start_time
    max_memory = usage_end.ru_maxrss - usage_start.ru_maxrss
    user_cpu_time = usage_end.ru_utime - usage_start.ru_utime
    system_cpu_time = usage_end.ru_stime - usage_start.ru_stime
    
    with open(sys.argv[1], 'w') as log_file:
        log_file.write(f"Elapsed time: {elapsed_time:.2f} seconds\n")
        log_file.write(f"Max memory usage: {max_memory} KB\n")
        log_file.write(f"User CPU time: {user_cpu_time:.2f} seconds\n")
        log_file.write(f"System CPU time: {system_cpu_time:.2f} seconds\n")
        log_file.write(f"Exit status: {result.returncode}\n")

if __name__ == "__main__":
    log_file = sys.argv[1]
    command = " ".join(sys.argv[2:])
    measure_resources(command)