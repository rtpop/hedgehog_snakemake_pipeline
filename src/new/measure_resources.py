import subprocess
import resource
import sys
import time
from contextlib import contextmanager


def _memory_unit():
    """Return the platform unit used by ru_maxrss."""
    return "KB" if sys.platform.startswith("linux") else "bytes"


def _write_measurements(log_file, start_time, usage_start, exit_status):
    usage_end = resource.getrusage(resource.RUSAGE_SELF)
    child_end = resource.getrusage(resource.RUSAGE_CHILDREN)
    elapsed_time = time.time() - start_time

    with open(log_file, "w") as output:
        output.write(f"Elapsed time: {elapsed_time:.2f} seconds\n")
        output.write(
            f"Max memory usage: {usage_end.ru_maxrss} {_memory_unit()}\n"
        )
        output.write(
            f"Max child memory usage: {child_end.ru_maxrss} {_memory_unit()}\n"
        )
        output.write(
            f"User CPU time: {usage_end.ru_utime - usage_start.ru_utime:.2f} seconds\n"
        )
        output.write(
            f"System CPU time: {usage_end.ru_stime - usage_start.ru_stime:.2f} seconds\n"
        )
        output.write(f"Exit status: {exit_status}\n")


@contextmanager
def measure_resources_context(log_file):
    """Measure work performed in the current Python process."""
    start_time = time.time()
    usage_start = resource.getrusage(resource.RUSAGE_SELF)
    try:
        yield
    except BaseException:
        _write_measurements(log_file, start_time, usage_start, 1)
        raise
    else:
        _write_measurements(log_file, start_time, usage_start, 0)

def measure_resources(command):
    start_time = time.time()
    usage_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    
    result = subprocess.run(command, shell=isinstance(command, str))
    
    usage_end = resource.getrusage(resource.RUSAGE_CHILDREN)
    end_time = time.time()
    
    elapsed_time = end_time - start_time
    max_memory = usage_end.ru_maxrss
    user_cpu_time = usage_end.ru_utime - usage_start.ru_utime
    system_cpu_time = usage_end.ru_stime - usage_start.ru_stime
    
    with open(sys.argv[1], 'w') as log_file:
        log_file.write(f"Elapsed time: {elapsed_time:.2f} seconds\n")
        log_file.write(f"Max memory usage: {max_memory} {_memory_unit()}\n")
        log_file.write(f"User CPU time: {user_cpu_time:.2f} seconds\n")
        log_file.write(f"System CPU time: {system_cpu_time:.2f} seconds\n")
        log_file.write(f"Exit status: {result.returncode}\n")

if __name__ == "__main__":
    log_file = sys.argv[1]
    command = sys.argv[2:]
    measure_resources(command)