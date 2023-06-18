cd "E:\__ONGOING\Projects\HORIZON\AstroSim\single body in gravitational well\"

g++  -o "one body.exe" "one body.cpp"

start /WAIT "running" "one body.exe"

start /WAIT python  "plot.py"

pause
