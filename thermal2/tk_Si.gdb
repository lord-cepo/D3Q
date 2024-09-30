file d3_tk.x
set cwd ../../Examples/Silicon
set args -in input.TK-sma
b sum_R2_sparse
r
p fc%dat(1)
