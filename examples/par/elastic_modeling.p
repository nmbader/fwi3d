# time parameters
courant=0.5
dt=-1
fmax=10
resampling=sinc
sinc_half_length=11
sub=0

# sources and receivers geometry
srcoord=./par/srcoord.txt

# source mechanism and receiver type
mt=0
fangle=1.5708
mxx=0
myy=0
mzz=0
mxy=0
mxz=0
myz=0
seismotype=1
gl=0

# boundary parameters
bc_top=1
bc_bottom=2
bc_left=2
bc_right=2
bc_front=2
bc_back=2
taper_top=0
taper_left=25
taper_right=25
taper_bottom=25
taper_front=25
taper_back=25
taper_strength=0.05

# model bounds
vpmin=0.2
vpmax=8
vsmin=0.1
vsmax=5
rhomin=0 
rhomax=8
deltamin=-0.5
deltamax=1
epsilonmin=-0.5
epsilonmax=1

# miscallenous
nthreads=24
verbose=3