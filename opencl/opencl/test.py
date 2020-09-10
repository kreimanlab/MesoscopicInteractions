import pyopencl as cl
from pyopencl import array
import numpy
import h5py
import time
 
if __name__ == "__main__":
    #vector = numpy.zeros((1, 4), cl.array.vec)
    #matrix = numpy.zeros((4, 4), cl.array.vec)
    vec_len = 3000
    mat_len = 3000
    mat_width = 20000
    n_compute = 2000
    vector = numpy.random.rand(1,vec_len).astype(numpy.float32)
    matrix = numpy.random.rand(mat_width,mat_len).astype(numpy.float32)
    for i in range(mat_width):
        matrix[i,:] += i + 1
     
    ## Step #1. Obtain an OpenCL platform.
    platform = cl.get_platforms()[0]
     
    ## It would be necessary to add some code to check the check the support for
    ## the necessary platform extensions with platform.extensions
     
    ## Step #2. Obtain a device id for at least one device (accelerator).
    device = platform.get_devices()[0]
     
    ## It would be necessary to add some code to check the check the support for
    ## the necessary device extensions with device.extensions
     
    ## Step #3. Create a context for the selected device.
    context = cl.Context([device])
     
    ## Step #4. Create the accelerator program from source code.
    ## Step #5. Build the program.
    ## Step #6. Create one or more kernels from the program functions.
    program = cl.Program(context, """
        __kernel void matrix_dot_vector(__global const float *matrix,
        __global const float *vector, __global float *result, const int vector_len)
        {
          int gid = get_global_id(0);

          float value = 0;
          for (unsigned int k = 0; k < vector_len; k++) {
              for (unsigned int l = 0; l < 10000; l++) {
                  value += acos(vector[k]) ;
              }
                  
              value += matrix[gid * vector_len + k] * vector[k];
          }
          result[gid] = value;
        }
        """).build()
     
    ## Step #7. Create a command queue for the target device.
    queue = cl.CommandQueue(context)
     
    ## Step #8. Allocate device memory and move input data from the host to the device memory.
    mem_flags = cl.mem_flags
    matrix_buf = cl.Buffer(context, mem_flags.READ_ONLY | mem_flags.COPY_HOST_PTR, hostbuf=matrix)
    vector_buf = cl.Buffer(context, mem_flags.READ_ONLY | mem_flags.COPY_HOST_PTR, hostbuf=vector)
    #vec_len_buf = cl.Buffer(context, mem_flags.READ_ONLY | mem_flags.COPY_HOST_PTR, hostbuf=numpy.array(vec_len,numpy.int))
    matrix_dot_vector = numpy.zeros(n_compute, numpy.float32) # mat_width
    destination_buf = cl.Buffer(context, mem_flags.WRITE_ONLY, matrix_dot_vector.nbytes)
     
    ## Step #9. Associate the arguments to the kernel with kernel object.
    ## Step #10. Deploy the kernel for device execution.
    program.matrix_dot_vector(queue, matrix_dot_vector.shape, None, matrix_buf, vector_buf, destination_buf, numpy.int32(vec_len))
     
    ## Step #11. Move the kernel’s output data to host memory.
    cl.enqueue_copy(queue, matrix_dot_vector, destination_buf)
     
    ## Step #12. Release context, program, kernels and memory.
    ## PyOpenCL performs this step for you, and therefore,
    ## you don't need to worry about cleanup code
     
    #print(matrix_dot_vector)
