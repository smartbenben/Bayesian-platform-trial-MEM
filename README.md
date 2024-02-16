main.R
apply different hybrid control methods (HY-MEM, HY-MAP, HY-rMAP) to simulated trials

main summary.R
calibrate each design assuming internal and external data is consistant by controling FWER at a target level, it produces the familiy-wise power and type I error rate of each
design under a varying level of internal-external data conflict in different scenarios

simulate.R
simulate internal data and external data

get.design.R
find the design parameters of each design under a prespecified FWER level

mem functions.R
implement the MEM approach and find a suitable choice for delta (find.mem.prior)
