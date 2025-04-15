## diamond-io

Implementation of [Diamond iO](https://eprint.iacr.org/2025/236)

### note

We currently have 2 different matrix implementation.
1. by calling `mmap()` syscall, use disk space as default storage
2. by using memory fully

#### Full Test (with test feature)

```
cargo test -r
```

#### Test iO (without test feature) 

this will remove helper logic + helper fields for test which is not require for iO security.


- dummy params
```
cargo test -r --test test_io_dummy_param --no-default-features (--features disk) -- --nocapture
```

- real params (by default ignored)
```
cargo test -r --test test_io_real_param --no-default-features (--features disk) -- --ignored --nocapture
```

- with memory profiler 
```
uv run memory_profile.py cargo test -r --test test_io_dummy_param (--features disk)  --no-default-features
```

#### Simulate Parameters
##### 1. Rust part: writing a rust program to simulate a norm corresponding to the circuit you want to obfuscate. You can refer to the `test_simulate_norm_final_bits_circuit` test function as an example. 

1.1. Fix the ring dimension `n`, the bits of each moduli of CRT `crt_bits`, the maximum number of moduli `crt_depth`, (implying that the bits of the derived modulus `q` is bounded by `crt_bits * crt_depth`,), and the base bits for decomposition `base_bits`.

1.2. Write a public circuit you want to obfuscate using `PolyCircuit`.

1.3. Construct `final_circuit` for your public circuit in the same manner as `test_simulate_norm_final_bits_circuit`.

1.4. Compute norms by calling the `final_circuit.simulate_bgg_norm` function. Note that the `unpacked_input_size` argument should be set to `1 + params.ring_dimension() * (the number of input polynomials to the public circuit except for the first 2 * log_base_q polynomials)`.

1.5. Store norms as a json file in the same manner as `test_simulate_norm_final_bits_circuit`.

##### 2. Python part: running `simulator/main.py` to get parameters satisfying correctness and security for your circuit.

2.1. Move the json file generated in Step 1.5 under the `simulator` directory.

2.2. Modify parameters after the line `if __name__ == "__main__":`. Specifically,
    - `secpar` denotes the minimum security parameter you want to guarantee.
    - `d` denotes the number of secret polynomials. We recommend setting it to `1` first and then increasing it when the simulator cannnot find any parameters with the current value.
    - `n`, `base_bits`, `crt_bits`, `max_crt_depth` should be the same as ones in Step 1.1.
    - `norms_path` is the file name you moved in Step 2.1.

2.3. Install sagemath if you have not installed it. Ref: https://doc.sagemath.org/html/en/installation/conda.html

2.4. Run `sage main.py`.

2.5. If the script is completed without any error, the found parameters are added to the last line in `simulator/params.log`. Among the parameters, `crt_depth` denotes the minimum number of moduli satisfying correctness and security, and `d`, `encoding_sigma`, `hardcoded_key_sigma`, `p_sigma`, and `switched_modulus` can be used for `ObfuscationParams`.

