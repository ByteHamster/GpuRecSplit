# GpuRecSplit / SimdRecSplit

We greatly improve the construction time of the [RecSplit](https://arxiv.org/abs/1910.06416) Minimal Perfect Hash Function using two orthogonal approaches.
*Rotation fitting* hashes the objects in each leaf to two sets and tries to combine them to a bijection by cyclically shifting one set to fill the holes in the other.
In addition, we harness parallelism on the level of bits, vectors, cores, and GPUs.
The code in this repository achieves significant speedups on SIMD machines and GPUs, compared
to the original [RecSplit implementation](https://github.com/vigna/sux/blob/master/sux/function/RecSplit.hpp).

[<img src="https://raw.githubusercontent.com/ByteHamster/GpuRecSplit/main/plots.png" alt="Plots preview">](https://arxiv.org/pdf/2212.09562)

| l | b | Method | Threads | B/Object | us/Object | Speedup |
|---:|---:|:---|---:|---:|---:|---:|
| 16 | 2000 | RecSplit [\[ALENEX'20\]](https://arxiv.org/abs/1910.06416) | 1 | 1.561 | 1152.6 |  |
| 16 | 2000 | SimdRecSplit | 1 | 1.561 | 139.2 | 8 |
| 16 | 2000 | SimdRecSplit | 16 | 1.562 | 28.1 | 40 |
| 16 | 2000 | GpuRecSplit | 4 | 1.562 | 1.5 | 763 |
| 18 | 50 | RecSplit [\[ALENEX'20\]](https://arxiv.org/abs/1910.06416) | 1 | 1.711 | 2919.5 |  |
| 18 | 50 | SimdRecSplit | 1 | 1.707 | 58.0 | 50 |
| 18 | 50 | SimdRecSplit | 16 | 1.709 | 12.1 | 241 |
| 18 | 50 | GpuRecSplit | 4 | 1.708 | 1.4 | 2072 |
| 24 | 2000 | GpuRecSplit | 4 | 1.524 | 633.9 |  |
    
### Library Usage

Clone (with submodules) this repo and add it to your `CMakeLists.txt`:

```
add_subdirectory(path/to/GpuRecSplit)
target_link_libraries(YourTarget PRIVATE RecSplit SIMDRecSplit GPURecSplit ShockHash) # or a subset of the targets
```

### Licensing
GpuRecSplit is licensed exactly like `libstdc++` (GPLv3 + GCC Runtime Library Exception), which essentially means you can use it everywhere, exactly like `libstdc++`.
You can find details in the [COPYING](/COPYING) and [COPYING.RUNTIME](/COPYING.RUNTIME) files.

If you use the project in an academic context or publication, please cite our paper:

```
@article{gpurecsplit2022,
  author    = {Dominik Bez and
        Florian Kurpicz and
        Hans{-}Peter Lehmann and
        Peter Sanders},
  title     = {High Performance Construction of
        RecSplit Based Minimal Perfect Hash Functions},
  journal   = {CoRR},
  volume    = {abs/2212.09562},
  year      = {2022},
  doi       = {10.48550/arXiv.2212.09562}
}
```
