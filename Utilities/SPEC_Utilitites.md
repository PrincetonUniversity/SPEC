# SPEC utilities

## Copy the `input` group from one file to another one using `h5py`:

```python
import h5py

f_src = h5py.File("demo_hdf5.sp.h5", "r")
f_tgt = h5py.File("demo_hdf5_2.sp.h5", "a")

# copy group '/input' from f_src into f_tgt
f_src.copy('/input', f_tgt)

f_tgt.close()
f_src.close()
```

