# GpuRecSplit / SimdRecSplit

Parallelization (Threads, SIMD, GPU) of the Minimal Perfect Hash Function
[RecSplit](https://arxiv.org/abs/1910.06416).
The code in this repository achieves significant speedups on SIMD machines and GPUs, compared
to the original [RecSplit implementation](https://github.com/vigna/sux/blob/master/sux/function/RecSplit.hpp).

<img src="https://raw.githubusercontent.com/ByteHamster/GpuRecSplit/main/plots.png" alt="Plots preview">

### Library Usage

Clone (with submodules) this repo and add it to your `CMakeLists.txt`:

```
add_subdirectory(path/to/GpuRecSplit)
target_link_libraries(YourTarget PRIVATE RecSplit SIMDRecSplit GPURecSplit ShockHash) # or a subset of the targets
```

### Licensing
GpuRecSplit is licensed exactly like `libstdc++` (GPLv3 + GCC Runtime Library Exception), which essentially means you can use it everywhere, exactly like `libstdc++`.
You can find details in the [COPYING](/COPYING) and [COPYING.RUNTIME](/COPYING.RUNTIME) files.
