## Add unregistered package
1. Open Julia REPL in project root directory.
2. Press `]` to enter package mode.
3. Use `activate .` to set environment.
4. Use `add https://github.com/yezhengkai/MINE_jll.jl` or `dev https://github.com/yezhengkai/MINE_jll.jl` to use `MINE_jll` package.

## Problems
1. Should we change the array's column-major memory to row-major memory?
2. Should we handle the memory of objects created by `ccall` or unsafe function?
