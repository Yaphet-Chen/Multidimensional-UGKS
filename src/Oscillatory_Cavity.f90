!--------------------------------------------------
!>Specify the KIND value
!--------------------------------------------------
module NumberKinds
    implicit none

    integer, parameter                                  :: KREAL = kind(0.d0)
    integer, parameter                                  :: KINT = kind(1)
end module NumberKinds

!--------------------------------------------------
!>Define the global constant variables
!--------------------------------------------------
module ConstantVariables
    use NumberKinds
    implicit none

    real(KREAL), parameter                              :: PI = 4.0*atan(1.0) !Pi
    real(KREAL), parameter                              :: SMV = tiny(0.0) !Small value to avoid 0/0
    real(KREAL), parameter                              :: UP = 1.0 !Used in sign() function
    
    !Reconstruction
    integer(KINT), parameter                            :: FIRST_ORDER = 0 !First order reconstruction
    integer(KINT), parameter                            :: LIMITER = 1 !VanleerLimiter reconstruction
    integer(KINT), parameter                            :: CENTRAL = 2 !Central difference reconstruction

    !Mesh type
    integer(KINT), parameter                            :: UNIFORM = 0 !Uniform mesh
    integer(KINT), parameter                            :: NONUNIFORM = 1 !Non-uniform mesh
    integer(KINT), parameter                            :: RAND = 2 !Random mesh

    !Output
    integer(KINT), parameter                            :: CENTER = 0 !Output solution as cell centered value
    integer(KINT), parameter                            :: POINTS = 1 !Output solution as point value

    !Quadrature method
    integer(KINT), parameter                            :: NEWTON = 0 !Newton–Cotes
    integer(KINT), parameter                            :: GAUSS = 1 !Gauss-Hermite
    integer(KINT), parameter                            :: TRAPEZOID = 2 !Trapezoidal rule

    !Boundary type
    integer(KINT), parameter                            :: KINETIC = 0 !Kinetic boundary condition
    integer(KINT), parameter                            :: MULTISCALE = 1 !Multiscale boundary condition

    !Flux type
    integer(KINT), parameter                            :: SPLITTING = 0 !Flux with direction splitting
    integer(KINT), parameter                            :: MULTIDIMENSION = 1 !Flux with multi-dimension
    
    !Direction
    integer(KINT), parameter                            :: IDIRC = 1 !I direction
    integer(KINT), parameter                            :: JDIRC = 2 !J direction

    !Rotation
    real(KREAL), parameter                              :: RN = 1 !No frame rotation
    real(KREAL), parameter                              :: RY = -1 !With frame rotation
end module ConstantVariables

!--------------------------------------------------
!>Define the mesh data structure
!--------------------------------------------------
module Mesh
    use NumberKinds
    implicit none

    !--------------------------------------------------
    !Basic derived type
    !--------------------------------------------------
    !Cell center
    type CellCenter
        !Geometry
        real(KREAL)                                     :: x,y !Cell center coordinates
        real(KREAL)                                     :: area !Cell area
        real(KREAL)                                     :: length(2) !Cell length in i and j direction
        !Flow field
        real(KREAL)                                     :: conVars(4) !Conservative variables at cell center: density, x-momentum, y-momentum, total energy
        real(KREAL), allocatable, dimension(:,:)        :: h,b !Distribution function
        real(KREAL), allocatable, dimension(:,:,:)      :: sh,sb !Slope of distribution function in i and j direction
    end type CellCenter

    !Cell interface
    type CellInterface
        !geometry
        real(KREAL)                                     :: length !Length of cell interface
        real(KREAL)                                     :: cosx,cosy !Directional cosine in global frame
        !Flow field
        real(KREAL)                                     :: conVars(4) !Conservative variables at cell interface: density, x-momentum, y-momentum, total energy
        !Flux
        real(KREAL)                                     :: flux(4) !Conservative variables flux at cell interface: density flux, x and y momentum flux, total energy flux
        real(KREAL), allocatable, dimension(:,:)        :: flux_h,flux_b !Flux of distribution function
    end type CellInterface

    !Grid geometry (node coordinates)
    type Grid
        real(KREAL)                                     :: x,y !Coordinates
    end type Grid

    !index method
    !---------------------------------
    !           (i,j+1)              |
    !      ----------------          |
    !      |              |          |
    !      |              |          |
    !      |              |          |
    ! (i,j)|     (i,j)    |(i+1,j)   |
    !      |      Cell    |          |
    !      |              |          |
    !      |              |          |
    !      ----------------          |
    !            (i,j)               |
    !---------------------------------
end module Mesh

module ControlParameters
    use ConstantVariables
    use Mesh
    implicit none

    !--------------------------------------------------
    !Variables to control the simulation
    !--------------------------------------------------
    integer(KINT), parameter                            :: RECONSTRUCTION_METHOD = LIMITER
    integer(KINT), parameter                            :: MESH_TYPE = NONUNIFORM
    integer(KINT), parameter                            :: QUADRATURE_TYPE = GAUSS
    integer(KINT), parameter                            :: OUTPUT_METHOD = CENTER
    integer(KINT), parameter                            :: BOUNDARY_TYPE = MULTISCALE
    integer(KINT), parameter                            :: FLUX_TYPE = MULTIDIMENSION
    real(KREAL), parameter                              :: CFL = 0.5 !CFL number
    integer(KINT), parameter                            :: MAX_ITER = 500000000 !Maximal iteration number
    real(KREAL), parameter                              :: EPS = 1.0E-7 !Convergence criteria
    real(KREAL)                                         :: simTime = 0.0 !Current simulation time
    integer(KINT)                                       :: iter = 1 !Number of iteration
    integer(KINT)                                       :: num_period = -1 !Number of accumulative oscillatory period
    real(KREAL)                                         :: dt !Global time step
    real(KREAL)                                         :: res(4) !Residual
    
    !Output control
    character(len=6), parameter                         :: HSTFILENAME = "Cavity" !History file name
    character(len=6), parameter                         :: RSTFILENAME = "Cavity" !Result file name
    character(len=6), parameter                         :: RESFILENAME = "Cavity" !Residual file name
    integer(KINT), parameter                            :: HSTFILE = 20 !History file ID
    integer(KINT), parameter                            :: RSTFILE = 21 !Result file ID
    integer(KINT), parameter                            :: RESFILE = 22 !Residual file ID

    !Gas propeties
    integer(KINT), parameter                            :: CK = 1 !Internal degree of freedom, here 1 denotes monatomic gas
    real(KREAL), parameter                              :: GAMMA = real(CK+4,KREAL)/real(CK+2,KREAL) !Ratio of specific heat
    real(KREAL), parameter                              :: OMEGA = 0.5 !Temperature dependence index in HS/VHS/VSS model
    real(KREAL), parameter                              :: PR = 2.0/3.0 !Prandtl number
    real(KREAL), parameter                              :: MA = 0.1 !Mach number
    real(KREAL), parameter                              :: ST = 2.0 !Strouhal number
    real(KREAL), parameter                              :: FRE = ST !Frequency omega
    real(KREAL), parameter                              :: TT = 2.0*PI/FRE !Peroid of Oscillation

    ! MU_REF determined by Kn number
    real(KREAL), parameter                              :: KN = 0.1 !Knudsen number in reference state
    real(KREAL), parameter                              :: ALPHA_REF = 1.0 !Coefficient in VHS model
    real(KREAL), parameter                              :: OMEGA_REF = 0.5 !Coefficient in VHS model
    real(KREAL), parameter                              :: MU_REF = 5.0*(ALPHA_REF+1.0)*(ALPHA_REF+2.0)*sqrt(PI)/(4.0*ALPHA_REF*(5.0-2.0*OMEGA_REF)*(7.0-2.0*OMEGA_REF))*KN !Viscosity coefficient in reference state

    ! MU_REF determined by Re number
    ! real(KREAL), parameter                              :: Re = 1000 !Reynolds number in reference state
    ! real(KREAL), parameter                              :: MU_REF = 0.15/Re !Viscosity coefficient in reference state

    !Geometry
    real(KREAL), parameter                              :: X_START = 0.0, X_END = 1.0, Y_START = 0.0, Y_END = 1.0 !Start point and end point in x, y direction 
    integer(KINT), parameter                            :: X_NUM = 31, Y_NUM = 31 !Points number in x, y direction
    integer(KINT), parameter                            :: IXMIN = 1 , IXMAX = X_NUM, IYMIN = 1 , IYMAX = Y_NUM !Cell index range
    integer(KINT), parameter                            :: N_GRID = (IXMAX-IXMIN+1)*(IYMAX-IYMIN+1) !Total number of cell
    
    !--------------------------------------------------
    !Discrete velocity space
    !--------------------------------------------------
    integer(KINT)                                       :: uNum = 16, vNum = 16 !Number of points in velocity space for u and v
    real(KREAL)                                         :: U_MIN = -5.0, U_MAX = +5.0, V_MIN = -5.0, V_MAX = +5.0 !Minimum and maximum micro velocity
    real(KREAL), allocatable, dimension(:,:)            :: uSpace,vSpace !Discrete velocity space for u and v
    real(KREAL), allocatable, dimension(:,:)            :: weight !Qudrature weight for discrete points in velocity space

    !--------------------------------------------------
    !Initial flow field
    !--------------------------------------------------
    !>Index method
    !---------------------------------
    !           (i,j+1)              |
    !      ----------------          |
    !      |              |          |
    !      |              |          |
    !      |              |          |
    ! (i,j)|     (i,j)    |(i+1,j)   |
    !      |      Cell    |          |
    !      |              |          |
    !      |              |          |
    !      ----------------          |
    !            (i,j)               |
    !---------------------------------
    type(CellCenter)                                    :: ctr(IXMIN:IXMAX,IYMIN:IYMAX) !Cell center
    type(CellInterface)                                 :: vface(IXMIN:IXMAX+1,IYMIN:IYMAX),hface(IXMIN:IXMAX,IYMIN:IYMAX+1) !Vertical and horizontal interfaces
    type(Grid)                                          :: geometry(IXMIN:IXMAX+1,IYMIN:IYMAX+1)

    !Initial condition (density, u-velocity, v-velocity, lambda=1/temperature)
    real(KREAL), parameter, dimension(4)                :: INIT_GAS = [1.0, 0.0, 0.0, 1.0]

    !Boundary condition (density, u-velocity, v-velocity, lambda=1/temperature)
    !------------------------------
    !              North
    !          ------------
    !          |           |
    !   West   |           |   East
    !          |           |
    !          ------------
    !              South
    !------------------------------
    real(KREAL), parameter                              :: U0 = MA*sqrt(GAMMA/2.0)
    real(KREAL), parameter, dimension(4)                :: BC_W = [1.0, 0.0, 0.0, 1.0] !West boundary
    real(KREAL), parameter, dimension(4)                :: BC_E = [1.0, 0.0, 0.0, 1.0] !East boundary
    real(KREAL), parameter, dimension(4)                :: BC_S = [1.0, 0.0, 0.0, 1.0] !South boundary
    real(KREAL), dimension(4)                           :: BC_N = [1.0, U0,  0.0, 1.0] !North boundary
end module ControlParameters

!--------------------------------------------------
!>Define some commonly used functions/subroutines
!--------------------------------------------------
module Tools
    use ControlParameters
    implicit none

contains
    !--------------------------------------------------
    !>Convert macro variables from global frame to local
    !>@param[in] w            :macro variables in global frame
    !>@param[in] cosx,cosy    :directional cosine
    !>@return    LocalFrame   :macro variables in local frame
    !--------------------------------------------------
    function LocalFrame(w,cosx,cosy)
        real(KREAL), intent(in)                         :: w(4)
        real(KREAL), intent(in)                         :: cosx,cosy
        real(KREAL)                                     :: LocalFrame(4)

        LocalFrame(1) = w(1)
        LocalFrame(2) = w(2)*cosx+w(3)*cosy
        LocalFrame(3) =-w(2)*cosy+w(3)*cosx
        LocalFrame(4) = w(4)
    end function LocalFrame

    !--------------------------------------------------
    !>Convert macro variables from local frame to global
    !>@param[in] w            :macro variables in local frame
    !>@param[in] cosx,cosy    :directional cosine
    !>@return    GlobalFrame  :macro variables in global frame
    !--------------------------------------------------
    function GlobalFrame(w,cosx,cosy)
        real(KREAL), intent(in)                         :: w(4)
        real(KREAL), intent(in)                         :: cosx,cosy
        real(KREAL)                                     :: GlobalFrame(4)

        GlobalFrame(1) = w(1)
        GlobalFrame(2) = w(2)*cosx-w(3)*cosy
        GlobalFrame(3) = w(2)*cosy+w(3)*cosx
        GlobalFrame(4) = w(4)
    end function GlobalFrame

    !--------------------------------------------------
    !>Convert primary variables to conservative variables in the equilibrium state
    !>@param[in] prim          :primary variables
    !>@return    GetConserved  :conservative variables
    !--------------------------------------------------
    function GetConserved(prim)
        real(KREAL), intent(in)                         :: prim(4) !Density, x-velocity, y-velocity, lambda=1/temperature
        real(KREAL)                                     :: GetConserved(4) !Density, x-momentum, y-momentum, total energy

        GetConserved(1) = prim(1)
        GetConserved(2) = prim(1)*prim(2)
        GetConserved(3) = prim(1)*prim(3)
        GetConserved(4) = 0.5*prim(1)/(prim(4)*(GAMMA-1.0))+0.5*prim(1)*(prim(2)**2+prim(3)**2)
    end function GetConserved

    !--------------------------------------------------
    !>Convert conservative variables to primary variables in the equilibrium state
    !>@param[in] w           :conservative variables
    !>@return    GetPrimary  :primary variables
    !--------------------------------------------------
    function GetPrimary(w)
        real(KREAL), intent(in)                         :: w(4) !Density, x-momentum, y-momentum, total energy
        real(KREAL)                                     :: GetPrimary(4) !Density, x-velocity, y-velocity, lambda=1/temperature

        GetPrimary(1) = w(1)
        GetPrimary(2) = w(2)/w(1)
        GetPrimary(3) = w(3)/w(1)
        GetPrimary(4) = 0.5*w(1)/((GAMMA-1.0)*(w(4)-0.5*(w(2)**2+w(3)**2)/w(1)))
    end function GetPrimary

    !--------------------------------------------------
    !>Obtain speed of sound
    !>@param[in] prim    :primary variables
    !>@return    GetSoundSpeed :speed of sound
    !--------------------------------------------------
    function GetSoundSpeed(prim)
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL)                                     :: GetSoundSpeed !Speed of sound

        GetSoundSpeed = sqrt(0.5*GAMMA/prim(4))
    end function GetSoundSpeed

    !--------------------------------------------------
    !>Get pressure
    !>@param[in] h,b          :distribution function
    !>@param[in] vn,vt        :normal and tangential velocity
    !>@param[in] prim         :primary variables
    !>@return    GetPressure  :pressure
    !--------------------------------------------------
    function GetPressure(h,b,vn,vt,prim)
    real(KREAL), dimension(:,:), intent(in)             :: h,b
    real(KREAL), dimension(:,:), intent(in)             :: vn,vt
    real(KREAL), intent(in)                             :: prim(4)
    real(KREAL)                                         :: GetPressure !pressure

    GetPressure = (sum(weight*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*b))/(CK+2)
    end function GetPressure

    !--------------------------------------------------
    !>Obtain discretized Maxwellian distribution
    !>@param[out] h,b   :distribution function
    !>@param[in]  vn,vt :normal and tangential velocity
    !>@param[in]  prim  :primary variables
    !--------------------------------------------------
    subroutine DiscreteMaxwell(h,b,vn,vt,prim)
        real(KREAL), dimension(:,:), intent(out)        :: h,b !Reduced distribution function
        real(KREAL), dimension(:,:), intent(in)         :: vn,vt !Normal and tangential velocity
        real(KREAL), intent(in)                         :: prim(4)

        h = prim(1)*(prim(4)/PI)*exp(-prim(4)*((vn-prim(2))**2+(vt-prim(3))**2))
        b = h*CK/(2.0*prim(4))
    end subroutine DiscreteMaxwell

    !--------------------------------------------------
    !>Calculate heat flux
    !>@param[in] h,b           :distribution function
    !>@param[in] vn,vt         :normal and tangential velocity
    !>@param[in] prim          :primary variables
    !>@return    GetHeatFlux   :heat flux in normal and tangential direction
    !--------------------------------------------------
    function GetHeatFlux(h,b,vn,vt,prim)
        real(KREAL), dimension(:,:), intent(in)         :: h,b
        real(KREAL), dimension(:,:), intent(in)         :: vn,vt
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL)                                     :: GetHeatFlux(2) !Heat flux in normal and tangential direction

        GetHeatFlux(1) = 0.5*(sum(weight*(vn-prim(2))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vn-prim(2))*b))
        GetHeatFlux(2) = 0.5*(sum(weight*(vt-prim(3))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vt-prim(3))*b))
    end function GetHeatFlux

    !--------------------------------------------------
    !>Calculate the Shakhov part H^+, B^+
    !>@param[in]  H,B           :Maxwellian distribution function
    !>@param[in]  vn,vt         :normal and tangential velocity
    !>@param[in]  qf            :heat flux
    !>@param[in]  prim          :primary variables
    !>@param[out] H_plus,B_plus :Shakhov part
    !--------------------------------------------------
    subroutine ShakhovPart(H,B,vn,vt,qf,prim,H_plus,B_plus)
        real(KREAL), dimension(:,:), intent(in)         :: H,B
        real(KREAL), dimension(:,:), intent(in)         :: vn,vt
        real(KREAL), intent(in)                         :: qf(2)
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL), dimension(:,:), intent(out)        :: H_plus,B_plus

        H_plus = 0.8*(1-PR)*prim(4)**2/prim(1)*&
                    ((vn-prim(2))*qf(1)+(vt-prim(3))*qf(2))*(2*prim(4)*((vn-prim(2))**2+(vt-prim(3))**2)+CK-5)*H
        B_plus = 0.8*(1-PR)*prim(4)**2/prim(1)*&
                    ((vn-prim(2))*qf(1)+(vt-prim(3))*qf(2))*(2*prim(4)*((vn-prim(2))**2+(vt-prim(3))**2)+CK-3)*B
    end subroutine ShakhovPart
    
    !--------------------------------------------------
    !>VanLeerLimiter for reconstruction of distrubution function
    !>@param[in]    leftCell  :the left cell
    !>@param[inout] midCell   :the middle cell
    !>@param[in]    rightCell :the right cell
    !>@param[in]    idx       :the index indicating i or j direction
    !--------------------------------------------------
    subroutine VanLeerLimiter(leftCell,midCell,rightCell,idx)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellCenter), intent(inout)                 :: midCell
        integer(KINT), intent(in)                       :: idx
        real(KREAL), allocatable, dimension(:,:)        :: sL,sR

        !Allocate array
        allocate(sL(uNum,vNum))
        allocate(sR(uNum,vNum))

        sL = (midCell%h-leftCell%h)/(0.5*(midCell%length(idx)+leftCell%length(idx)))
        sR = (rightCell%h-midCell%h)/(0.5*(rightCell%length(idx)+midCell%length(idx)))
        midCell%sh(:,:,idx) = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

        sL = (midCell%b-leftCell%b)/(0.5*(midCell%length(idx)+leftCell%length(idx)))
        sR = (rightCell%b-midCell%b)/(0.5*(rightCell%length(idx)+midCell%length(idx)))
        midCell%sb(:,:,idx) = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

        !Deallocate array
        deallocate(sL)
        deallocate(sR)
    end subroutine VanLeerLimiter
end module Tools

!--------------------------------------------------
!>Flux calculation
!--------------------------------------------------
module Flux
    use Tools
    implicit none
    interface CalcFlux
        module procedure CalcFlux_MultiDimension,CalcFlux_DirectionSplitting
    end interface
    integer(KREAL), parameter                           :: MNUM = 6 !Number of normal velocity moments
    integer(KREAL), parameter                           :: MTUM = 5 !Number of tangential velocity moments

contains
    !--------------------------------------------------
    !>Calculate conVars at cell interface in local frame
    !>@param[in]    leftCell  :cell left to the target interface
    !>@param[inout] face      :the target interface
    !>@param[in]    rightCell :cell right to the target interface
    !>@param[in]    idx       :index indicating i or j direction
    !--------------------------------------------------
    subroutine CalcFaceConvars(leftCell,face,rightCell)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellInterface), intent(inout)              :: face

        face%conVars = 0.5*LocalFrame(rightCell%conVars+leftCell%conVars,face%cosx,face%cosy)
    end subroutine CalcFaceConvars

    !--------------------------------------------------
    !>Calculate flux of inner interface with multi-dimension
    !>@param[in]    leftCell  :cell left to the target interface
    !>@param[inout] face      :the target interface
    !>@param[in]    rightCell :cell right to the target interface
    !>@param[in]    idx       :index indicating i or j direction
    !--------------------------------------------------
    subroutine CalcFlux_MultiDimension(leftCell,face,rightCell,idx,face_U,face_D,idb)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellInterface), intent(inout)              :: face
        type(CellInterface), intent(in)                 :: face_U,face_D !Upper and Lower interface
        integer(KINT), intent(in)                       :: idx,idb !indicator for direction and boundary
        real(KREAL), allocatable, dimension(:,:)        :: vn,vt !normal and tangential micro velocity
        real(KREAL), allocatable, dimension(:,:)        :: h,b !Distribution function at the interface
        real(KREAL), allocatable, dimension(:,:)        :: H0,B0 !Maxwellian distribution function
        real(KREAL), allocatable, dimension(:,:)        :: H_plus,B_plus !Shakhov part of the equilibrium distribution
        real(KREAL), allocatable, dimension(:,:)        :: shn,sbn !Slope of distribution function at the interface
        real(KREAL), allocatable, dimension(:,:)        :: sht,sbt !Tangential slope of distribution function at the interface
        real(KREAL), allocatable, dimension(:,:)        :: delta !Heaviside step function
        real(KREAL)                                     :: prim(4) !Primary variables at the interface
        real(KREAL)                                     :: qf(2) !Heat flux in normal and tangential direction
        real(KREAL)                                     :: sw_n(4),sw_t(4) !Slope of conVars at normal and tangential direction
        real(KREAL)                                     :: a_slope(4),b_slope(4),aT(4) !Micro slope of Maxwellian distribution at normal tangential and time.
        real(KREAL)                                     :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM),Mv(0:MTUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<v^m>,<\xi^l>
        real(KREAL)                                     :: Mau0(4),Mau(4),Mbv(4),MauT(4) !<u\psi>,<a*u^n*\psi>,<b*u*v*\psi>,<A*u*\psi>
        real(KREAL)                                     :: tau !Collision time
        real(KREAL)                                     :: Mt(5) !Some time integration terms

        !--------------------------------------------------
        !Prepare
        !--------------------------------------------------
        !Allocate array
        allocate(vn(uNum,vNum))
        allocate(vt(uNum,vNum))
        allocate(delta(uNum,vNum))
        allocate(h(uNum,vNum))
        allocate(b(uNum,vNum))
        allocate(shn(uNum,vNum))
        allocate(sbn(uNum,vNum))
        allocate(sht(uNum,vNum))
        allocate(sbt(uNum,vNum))
        allocate(H0(uNum,vNum))
        allocate(B0(uNum,vNum))
        allocate(H_plus(uNum,vNum))
        allocate(B_plus(uNum,vNum))

        !Convert the velocity space to local frame
        vn = uSpace*face%cosx+vSpace*face%cosy
        vt =-uSpace*face%cosy+vSpace*face%cosx

        !Heaviside step function
        delta = (sign(UP,vn)+1.0)/2.0

        !--------------------------------------------------
        !Reconstruct initial distribution at interface
        !--------------------------------------------------
        h = (leftCell%h+0.5*leftCell%length(idx)*leftCell%sh(:,:,idx))*delta+&
            (rightCell%h-0.5*rightCell%length(idx)*rightCell%sh(:,:,idx))*(1-delta)
        b = (leftCell%b+0.5*leftCell%length(idx)*leftCell%sb(:,:,idx))*delta+&
            (rightCell%b-0.5*rightCell%length(idx)*rightCell%sb(:,:,idx))*(1-delta)
        shn = leftCell%sh(:,:,idx)*delta+rightCell%sh(:,:,idx)*(1-delta)
        sbn = leftCell%sb(:,:,idx)*delta+rightCell%sb(:,:,idx)*(1-delta)

        !Tangential part
        if (idx==1) then
            sht = leftCell%sh(:,:,2)*delta+rightCell%sh(:,:,2)*(1-delta)
            sbt = leftCell%sb(:,:,2)*delta+rightCell%sb(:,:,2)*(1-delta)
        else
            sht = -(leftCell%sh(:,:,1)*delta+rightCell%sh(:,:,1)*(1-delta))
            sbt = -(leftCell%sb(:,:,1)*delta+rightCell%sb(:,:,1)*(1-delta))
        end if
        !--------------------------------------------------
        !Obtain macroscopic variables
        !--------------------------------------------------
        !Obtain primary variables at interface
        prim = GetPrimary(face%conVars) !face%conVars is in local frame already

        !--------------------------------------------------
        !Calculate a_slope
        !--------------------------------------------------
        sw_n = LocalFrame(rightCell%conVars-leftCell%conVars,face%cosx,face%cosy)/(0.5*leftCell%length(idx)+0.5*rightCell%length(idx)) !normal slope of face%conVars
        a_slope = MicroSlope(prim,sw_n) !Calculate a^L

        sw_t = (face_U%conVars-face_D%conVars)/(0.5*face_U%length+idb*face%length+0.5*face_D%length) !right slope of face%conVars
        b_slope = MicroSlope(prim,sw_t) !Calculate b_slope

        !--------------------------------------------------
        !Calculate time slope of conVars and A
        !--------------------------------------------------
        !<u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
        call CalcMoment(prim,Mu,Mv,Mxi,MuL,MuR) 

        Mau = Moment_auvxi(a_slope,Mu,Mv,Mxi,1,0) !<aL*u*\psi>
        Mbv = Moment_auvxi(b_slope,Mu,Mv,Mxi,0,1) !<b*u*v*\psi>

        aT = MicroSlope(prim,-prim(1)*(Mau+Mbv)) !Calculate A

        !--------------------------------------------------
        !Calculate collision time and some time integration terms
        !--------------------------------------------------
        tau = GetTau(prim)

        Mt(4) = tau*(1.0-exp(-dt/tau))
        Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
        Mt(1) = dt-Mt(4)
        Mt(2) = -tau*Mt(1)+Mt(5) 
        Mt(3) = 0.5*dt**2-tau*Mt(1)

        !--------------------------------------------------
        !Calculate the flux of conservative variables related to g0
        !--------------------------------------------------
        Mau0 = Moment_uvxi(Mu,Mv,Mxi,1,0,0) !<u*\psi>
        Mau = Moment_auvxi(a_slope,Mu,Mv,Mxi,2,0) !<a*u^2*\psi>
        Mbv = Moment_auvxi(b_slope,Mu,Mv,Mxi,1,1) !<b*u*v*\psi>
        MauT = Moment_auvxi(aT,Mu,Mv,Mxi,1,0) !<A*u*\psi>

        face%flux = Mt(1)*prim(1)*Mau0+Mt(2)*prim(1)*(Mau+Mbv)+Mt(3)*prim(1)*MauT

        !--------------------------------------------------
        !Calculate the flux of conservative variables related to g+ and f0
        !--------------------------------------------------
        !Maxwellian distribution H0 and B0
        call DiscreteMaxwell(H0,B0,vn,vt,prim)
    
        !Calculate heat flux
        qf = GetHeatFlux(h,b,vn,vt,prim)

        !Shakhov part H+ and B+
        call ShakhovPart(H0,B0,vn,vt,qf,prim,H_plus,B_plus)

        !Conservative flux related to g+ and f0
        face%flux(1) = face%flux(1)+Mt(1)*sum(weight*vn*H_plus)+Mt(4)*sum(weight*vn*h)-Mt(5)*sum(weight*vn*(vn*shn+vt*sht))
        face%flux(2) = face%flux(2)+Mt(1)*sum(weight*vn*vn*H_plus)+Mt(4)*sum(weight*vn*vn*h)-Mt(5)*sum(weight*vn**2*(vn*shn+vt*sht))
        face%flux(3) = face%flux(3)+Mt(1)*sum(weight*vt*vn*H_plus)+Mt(4)*sum(weight*vt*vn*h)-Mt(5)*sum(weight*vt*vn*(vn*shn+vt*sht))
        face%flux(4) = face%flux(4)+&
                        Mt(1)*0.5*sum(weight*vn*((vn**2+vt**2)*H_plus+B_plus))+&
                        Mt(4)*0.5*sum(weight*vn*((vn**2+vt**2)*h+b))-&
                        Mt(5)*0.5*sum(weight*vn*((vn**2+vt**2)*(vn*shn+vt*sht)+(vn*sbn+vt*sht)))

        !--------------------------------------------------
        !Calculate flux of distribution function
        !--------------------------------------------------
        face%flux_h = Mt(1)*vn*(H0+H_plus)+&
                        Mt(2)*vn**2*(a_slope(1)*H0+a_slope(2)*vn*H0+a_slope(3)*vt*H0+0.5*a_slope(4)*((vn**2+vt**2)*H0+B0))+&
                        Mt(2)*vn*vt*(b_slope(1)*H0+b_slope(2)*vn*H0+b_slope(3)*vt*H0+0.5*b_slope(4)*((vn**2+vt**2)*H0+B0))+&
                        Mt(3)*vn*(aT(1)*H0+aT(2)*vn*H0+aT(3)*vt*H0+0.5*aT(4)*((vn**2+vt**2)*H0+B0))+&
                        Mt(4)*vn*h-Mt(5)*vn*(vn*shn+vt*sht)

        face%flux_b = Mt(1)*vn*(B0+B_plus)+&
                        Mt(2)*vn**2*(a_slope(1)*B0+a_slope(2)*vn*B0+a_slope(3)*vt*B0+0.5*a_slope(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                        Mt(2)*vn*vt*(b_slope(1)*B0+b_slope(2)*vn*B0+b_slope(3)*vt*B0+0.5*b_slope(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                        Mt(3)*vn*(aT(1)*B0+aT(2)*vn*B0+aT(3)*vt*B0+0.5*aT(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                        Mt(4)*vn*b-Mt(5)*vn*(vn*sbn+vt*sht)

        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !Convert to global frame
        face%flux = GlobalFrame(face%flux,face%cosx,face%cosy) 
        !Total flux
        face%flux = face%length*face%flux
        face%flux_h = face%length*face%flux_h
        face%flux_b = face%length*face%flux_b

        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Deallocate array
        deallocate(vn)
        deallocate(vt)
        deallocate(delta)
        deallocate(h)
        deallocate(b)
        deallocate(shn)
        deallocate(sbn)
        deallocate(sht)
        deallocate(sbt)
        deallocate(H0)
        deallocate(B0)
        deallocate(H_plus)
        deallocate(B_plus)
    end subroutine CalcFlux_MultiDimension
    
    !--------------------------------------------------
    !>Calculate flux of inner interface with direction splitting
    !>@param[in]    leftCell  :cell left to the target interface
    !>@param[inout] face      :the target interface
    !>@param[in]    rightCell :cell right to the target interface
    !>@param[in]    idx       :index indicating i or j direction
    !--------------------------------------------------
    subroutine CalcFlux_DirectionSplitting(leftCell,face,rightCell,idx)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellInterface), intent(inout)              :: face
        integer(KINT), intent(in)                       :: idx
        real(KREAL), allocatable, dimension(:,:)        :: vn,vt !normal and tangential micro velocity
        real(KREAL), allocatable, dimension(:,:)        :: h,b !Distribution function at the interface
        real(KREAL), allocatable, dimension(:,:)        :: H0,B0 !Maxwellian distribution function
        real(KREAL), allocatable, dimension(:,:)        :: H_plus,B_plus !Shakhov part of the equilibrium distribution
        real(KREAL), allocatable, dimension(:,:)        :: sh,sb !Slope of distribution function at the interface
        integer(KINT), allocatable, dimension(:,:)      :: delta !Heaviside step function
        real(KREAL)                                     :: conVars(4),prim(4) !Conservative and primary variables at the interface
        real(KREAL)                                     :: qf(2) !Heat flux in normal and tangential direction
        real(KREAL)                                     :: sw(4) !Slope of conVars
        real(KREAL)                                     :: aL(4),aR(4),aT(4) !Micro slope of Maxwellian distribution, left,right and time.
        real(KREAL)                                     :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM),Mv(0:MTUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<v^m>,<\xi^l>
        real(KREAL)                                     :: Mau0(4),MauL(4),MauR(4),MauT(4) !<u\psi>,<aL*u^n*\psi>,<aR*u^n*\psi>,<A*u*\psi>
        real(KREAL)                                     :: tau !Collision time
        real(KREAL)                                     :: Mt(5) !Some time integration terms
        integer(KINT)                                   :: i,j

        !--------------------------------------------------
        !Prepare
        !--------------------------------------------------
        !Allocate array
        allocate(vn(uNum,vNum))
        allocate(vt(uNum,vNum))
        allocate(delta(uNum,vNum))
        allocate(h(uNum,vNum))
        allocate(b(uNum,vNum))
        allocate(sh(uNum,vNum))
        allocate(sb(uNum,vNum))
        allocate(H0(uNum,vNum))
        allocate(B0(uNum,vNum))
        allocate(H_plus(uNum,vNum))
        allocate(B_plus(uNum,vNum))

        !Convert the velocity space to local frame
        vn = uSpace*face%cosx+vSpace*face%cosy
        vt =-uSpace*face%cosy+vSpace*face%cosx

        !Heaviside step function
        delta = (sign(UP,vn)+1)/2

        !--------------------------------------------------
        !Reconstruct initial distribution at interface
        !--------------------------------------------------
        h = (leftCell%h+0.5*leftCell%length(idx)*leftCell%sh(:,:,idx))*delta+&
            (rightCell%h-0.5*rightCell%length(idx)*rightCell%sh(:,:,idx))*(1-delta)
        b = (leftCell%b+0.5*leftCell%length(idx)*leftCell%sb(:,:,idx))*delta+&
            (rightCell%b-0.5*rightCell%length(idx)*rightCell%sb(:,:,idx))*(1-delta)
        sh = leftCell%sh(:,:,idx)*delta+rightCell%sh(:,:,idx)*(1-delta)
        sb = leftCell%sb(:,:,idx)*delta+rightCell%sb(:,:,idx)*(1-delta)

        !--------------------------------------------------
        !Obtain macroscopic variables
        !--------------------------------------------------
        !Conservative variables conVars at interface
        conVars(1) = sum(weight*h)
        conVars(2) = sum(weight*vn*h)
        conVars(3) = sum(weight*vt*h)
        conVars(4) = 0.5*(sum(weight*(vn**2+vt**2)*h)+sum(weight*b))

        !Convert to primary variables
        prim = GetPrimary(conVars)

        !--------------------------------------------------
        !Calculate a^L,a^R
        !--------------------------------------------------
        sw = (conVars-LocalFrame(leftCell%conVars,face%cosx,face%cosy))/(0.5*leftCell%length(idx)) !left slope of conVars
        aL = MicroSlope(prim,sw) !Calculate a^L

        sw = (LocalFrame(rightCell%conVars,face%cosx,face%cosy)-conVars)/(0.5*rightCell%length(idx)) !right slope of conVars
        aR = MicroSlope(prim,sw) !Calculate a^R

        !--------------------------------------------------
        !Calculate time slope of conVars and A
        !--------------------------------------------------
        !<u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
        call CalcMoment(prim,Mu,Mv,Mxi,MuL,MuR) 

        MauL = Moment_auvxi(aL,MuL,Mv,Mxi,1,0) !<aL*u*\psi>_{>0}
        MauR = Moment_auvxi(aR,MuR,Mv,Mxi,1,0) !<aR*u*\psi>_{<0}

        sw = -prim(1)*(MauL+MauR) !Time slope of conVars
        aT = MicroSlope(prim,sw) !Calculate A

        !--------------------------------------------------
        !Calculate collision time and some time integration terms
        !--------------------------------------------------
        tau = GetTau(prim)

        Mt(4) = tau*(1.0-exp(-dt/tau))
        Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
        Mt(1) = dt-Mt(4)
        Mt(2) = -tau*Mt(1)+Mt(5) 
        Mt(3) = 0.5*dt**2-tau*Mt(1)

        !--------------------------------------------------
        !Calculate the flux of conservative variables related to g0
        !--------------------------------------------------
        Mau0 = Moment_uvxi(Mu,Mv,Mxi,1,0,0) !<u*\psi>
        MauL = Moment_auvxi(aL,MuL,Mv,Mxi,2,0) !<aL*u^2*\psi>_{>0}
        MauR = Moment_auvxi(aR,MuR,Mv,Mxi,2,0) !<aR*u^2*\psi>_{<0}
        MauT = Moment_auvxi(aT,Mu,Mv,Mxi,1,0) !<A*u*\psi>

        face%flux = Mt(1)*prim(1)*Mau0+Mt(2)*prim(1)*(MauL+MauR)+Mt(3)*prim(1)*MauT

        !--------------------------------------------------
        !Calculate the flux of conservative variables related to g+ and f0
        !--------------------------------------------------
        !Maxwellian distribution H0 and B0
        call DiscreteMaxwell(H0,B0,vn,vt,prim)
    
        !Calculate heat flux
        qf = GetHeatFlux(h,b,vn,vt,prim) 

        !Shakhov part H+ and B+
        call ShakhovPart(H0,B0,vn,vt,qf,prim,H_plus,B_plus)

        !Conservative flux related to g+ and f0
        face%flux(1) = face%flux(1)+Mt(1)*sum(weight*vn*H_plus)+Mt(4)*sum(weight*vn*h)-Mt(5)*sum(weight*vn**2*sh)
        face%flux(2) = face%flux(2)+Mt(1)*sum(weight*vn*vn*H_plus)+Mt(4)*sum(weight*vn*vn*h)-Mt(5)*sum(weight*vn*vn**2*sh)
        face%flux(3) = face%flux(3)+Mt(1)*sum(weight*vt*vn*H_plus)+Mt(4)*sum(weight*vt*vn*h)-Mt(5)*sum(weight*vt*vn**2*sh)
        face%flux(4) = face%flux(4)+&
                        Mt(1)*0.5*(sum(weight*vn*(vn**2+vt**2)*H_plus)+sum(weight*vn*B_plus))+&
                        Mt(4)*0.5*(sum(weight*vn*(vn**2+vt**2)*h)+sum(weight*vn*b))-&
                        Mt(5)*0.5*(sum(weight*vn**2*(vn**2+vt**2)*sh)+sum(weight*vn**2*sb))

        !--------------------------------------------------
        !Calculate flux of distribution function
        !--------------------------------------------------
        face%flux_h = Mt(1)*vn*(H0+H_plus)+&
                        Mt(2)*vn**2*(aL(1)*H0+aL(2)*vn*H0+aL(3)*vt*H0+0.5*aL(4)*((vn**2+vt**2)*H0+B0))*delta+&
                        Mt(2)*vn**2*(aR(1)*H0+aR(2)*vn*H0+aR(3)*vt*H0+0.5*aR(4)*((vn**2+vt**2)*H0+B0))*(1-delta)+&
                        Mt(3)*vn*(aT(1)*H0+aT(2)*vn*H0+aT(3)*vt*H0+0.5*aT(4)*((vn**2+vt**2)*H0+B0))+&
                        Mt(4)*vn*h-Mt(5)*vn**2*sh

        face%flux_b = Mt(1)*vn*(B0+B_plus)+&
                        Mt(2)*vn**2*(aL(1)*B0+aL(2)*vn*B0+aL(3)*vt*B0+0.5*aL(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*delta+&
                        Mt(2)*vn**2*(aR(1)*B0+aR(2)*vn*B0+aR(3)*vt*B0+0.5*aR(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))*(1-delta)+&
                        Mt(3)*vn*(aT(1)*B0+aT(2)*vn*B0+aT(3)*vt*B0+0.5*aT(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
                        Mt(4)*vn*b-Mt(5)*vn**2*sb
        
        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !Convert to global frame
        face%flux = GlobalFrame(face%flux,face%cosx,face%cosy) 
        !Total flux
        face%flux = face%length*face%flux
        face%flux_h = face%length*face%flux_h
        face%flux_b = face%length*face%flux_b

        
        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Deallocate array
        deallocate(vn)
        deallocate(vt)
        deallocate(delta)
        deallocate(h)
        deallocate(b)
        deallocate(sh)
        deallocate(sb)
        deallocate(H0)
        deallocate(B0)
        deallocate(H_plus)
        deallocate(B_plus)
    end subroutine CalcFlux_DirectionSplitting

    !--------------------------------------------------
    !>Calculate flux of boundary interface, assuming left wall
    !>@param[in]    bc   :boundary condition
    !>@param[inout] face :the boundary interface
    !>@param[in]    cell :cell next to the boundary interface
    !>@param[in]    idx  :index indicating i or j direction
    !>@param[in]    rot  :indicating rotation
    !--------------------------------------------------
    subroutine CalcFluxBoundary(bc,face,cell,idx,rot)
        real(KREAL), intent(in)                         :: bc(4) !Primary variables at boundary
        type(CellInterface), intent(inout)              :: face
        type(CellCenter), intent(in)                    :: cell
        integer(KINT), intent(in)                       :: idx
        real(KREAL), intent(in)                         :: rot

        if (BOUNDARY_TYPE==KINETIC) then
            call KineticFluxBoundary(bc,face,cell,idx,rot)
        elseif (BOUNDARY_TYPE==MULTISCALE) then
            call MultiscaleFluxBoundary(bc,face,cell,idx,rot)
        else
            stop "Error in BOUNDARY_TYPE!"
        end if
    end subroutine CalcFluxBoundary

    !--------------------------------------------------
    !>Calculate kinetic flux of boundary interface, assuming left wall
    !>@param[in]    bc   :boundary condition
    !>@param[inout] face :the boundary interface
    !>@param[in]    cell :cell next to the boundary interface
    !>@param[in]    idx  :index indicating i or j direction
    !>@param[in]    rot  :indicating rotation
    !--------------------------------------------------
    subroutine KineticFluxBoundary(bc,face,cell,idx,rot)
        real(KREAL), intent(in)                         :: bc(4) !Primary variables at boundary
        type(CellInterface), intent(inout)              :: face
        type(CellCenter), intent(in)                    :: cell
        integer(KINT), intent(in)                       :: idx
        real(KREAL), intent(in)                         :: rot
        real(KREAL), allocatable, dimension(:,:)        :: vn,vt !Normal and tangential micro velocity
        real(KREAL), allocatable, dimension(:,:)        :: h,b !Reduced distribution function
        real(KREAL), allocatable, dimension(:,:)        :: H0,B0 !Maxwellian distribution function at the wall
        real(KREAL), allocatable, dimension(:,:)        :: delta !Heaviside step function
        real(KREAL)                                     :: prim(4) !boundary condition in local frame
        real(KREAL)                                     :: incidence,reflection
        !--------------------------------------------------
        !prepare
        !--------------------------------------------------
        !allocate array
        allocate(vn(uNum,vNum))
        allocate(vt(uNum,vNum))
        allocate(delta(uNum,vNum))
        allocate(h(uNum,vNum))
        allocate(b(uNum,vNum))
        allocate(H0(uNum,vNum))
        allocate(B0(uNum,vNum))

        !Convert the micro velocity to local frame
        vn = uSpace*face%cosx+vSpace*face%cosy
        vt =-uSpace*face%cosy+vSpace*face%cosx

        !Heaviside step function. The rotation accounts for the right wall
        delta = (sign(UP,vn)*rot+1.0)/2.0

        !Boundary condition in local frame
        prim = LocalFrame(bc,face%cosx,face%cosy)

        !--------------------------------------------------
        !Obtain h^{in} and b^{in}, rotation accounts for the right wall
        !--------------------------------------------------
        h = cell%h-rot*0.5*cell%length(idx)*cell%sh(:,:,idx)
        b = cell%b-rot*0.5*cell%length(idx)*cell%sb(:,:,idx)
        
        !--------------------------------------------------
        !Calculate wall density and Maxwellian distribution
        !--------------------------------------------------
        incidence = sum(weight*vn*h*(1-delta))
        reflection = (prim(4)/PI)*sum(weight*vn*exp(-prim(4)*((vn-prim(2))**2+(vt-prim(3))**2))*delta)

        prim(1) = -incidence/reflection

        call DiscreteMaxwell(H0,B0,vn,vt,prim)
        
        !--------------------------------------------------
        !Distribution function at the boundary interface
        !--------------------------------------------------
        h = H0*delta+h*(1-delta)
        b = B0*delta+b*(1-delta)
        
        !--------------------------------------------------
        !Calculate flux
        !--------------------------------------------------
        face%flux(1) = sum(weight*vn*h)
        face%flux(2) = sum(weight*vn*vn*h)
        face%flux(3) = sum(weight*vn*vt*h)
        face%flux(4) = 0.5*sum(weight*vn*((vn**2+vt**2)*h+b))

        face%flux_h = vn*h
        face%flux_b = vn*b

        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !Convert to global frame
        face%flux = GlobalFrame(face%flux,face%cosx,face%cosy)

        !Total flux
        face%flux = dt*face%length*face%flux
        face%flux_h = dt*face%length*face%flux_h
        face%flux_b = dt*face%length*face%flux_b
        
        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Deallocate array
        deallocate(vn)
        deallocate(vt)
        deallocate(delta)
        deallocate(h)
        deallocate(b)
        deallocate(H0)
        deallocate(B0)
    end subroutine KineticFluxBoundary

    !--------------------------------------------------
    !>Calculate kinetic flux of boundary interface, assuming left wall
    !>@param[in]    bc   :boundary condition
    !>@param[inout] face :the boundary interface
    !>@param[in]    cell :cell next to the boundary interface
    !>@param[in]    idx  :index indicating i or j direction
    !>@param[in]    rot  :indicating rotation
    !--------------------------------------------------
    subroutine MultiscaleFluxBoundary(bc,face,cell,idx,rot)
        real(KREAL), intent(in)                         :: bc(4) !Primary variables at boundary
        type(CellInterface), intent(inout)              :: face
        type(CellCenter), intent(in)                    :: cell
        integer(KINT), intent(in)                       :: idx
        real(KREAL), intent(in)                         :: rot
        real(KREAL), allocatable, dimension(:,:)        :: vn,vt !Normal and tangential micro velocity
        real(KREAL), allocatable, dimension(:,:)        :: h,b !Reduced non-equlibrium distribution function at interface
        real(KREAL), allocatable, dimension(:,:)        :: H_g,B_g !Maxwellian distribution function g_{i+1/2}
        real(KREAL), allocatable, dimension(:,:)        :: H_w,B_w !Maxwellian distribution function at the wall
        real(KREAL), allocatable, dimension(:,:)        :: delta !Heaviside step function
        real(KREAL)                                     :: prim(4)
        real(KREAL)                                     :: prim_g(4) !Primary variables of g_{i+1/2}
        real(KREAL)                                     :: prim_w(4) !boundary condition in local frame
        real(KREAL)                                     :: incidence,reflection
        real(KREAL)                                     :: tau,T1,T4,T5 !Some time integration coefficients

        !--------------------------------------------------
        !prepare
        !--------------------------------------------------
        !allocate array
        allocate(vn(uNum,vNum))
        allocate(vt(uNum,vNum))
        allocate(delta(uNum,vNum))
        allocate(h(uNum,vNum))
        allocate(b(uNum,vNum))
        allocate(H_g(uNum,vNum))
        allocate(B_g(uNum,vNum))
        allocate(H_w(uNum,vNum))
        allocate(B_w(uNum,vNum))

        !Convert the micro velocity to local frame
        vn = uSpace*face%cosx+vSpace*face%cosy
        vt =-uSpace*face%cosy+vSpace*face%cosx

        !Heaviside step function. The rotation accounts for the right wall
        delta = (sign(UP,vn)*rot+1.0)/2.0

        !--------------------------------------------------
        !Reconstruct non-equlibrium distribution at interface
        !--------------------------------------------------
        h = cell%h-rot*0.5*cell%length(idx)*cell%sh(:,:,idx)
        b = cell%b-rot*0.5*cell%length(idx)*cell%sb(:,:,idx)

        !Calculate primary variable of g_{i+1/2} in local frame
        prim_w = LocalFrame(bc,face%cosx,face%cosy)
        prim_g = prim_w
        prim = LocalFrame(GetPrimary(cell%conVars),face%cosx,face%cosy)
        prim_g(1) = prim(1)/prim(4)*prim_w(4) !keep consistent pressure inside cell

        !--------------------------------------------------
        !Calculate collision time and some time integration terms
        !--------------------------------------------------
        tau = GetTau(prim_g)

        T4 = tau*(1.0-exp(-dt/tau))
        T5 = -tau*dt*exp(-dt/tau)+tau*T4
        T1 = dt-T4

        !--------------------------------------------------
        !Calculate wall density and Maxwellian distribution
        !--------------------------------------------------
        call DiscreteMaxwell(H_g,B_g,vn,vt,prim_g)

        incidence = T1*sum(weight*vn*H_g*(1-delta))&
                    +T4*sum(weight*vn*h*(1-delta))-T5*sum(weight*vn*vn*cell%sh(:,:,idx)*(1-delta))
        reflection = dt*(prim_w(4)/PI)*sum(weight*vn*exp(-prim_w(4)*((vn-prim_w(2))**2+(vt-prim_w(3))**2))*delta)
        prim_w(1) = -incidence/reflection

        call DiscreteMaxwell(H_w,B_w,vn,vt,prim_w)

        !--------------------------------------------------
        !Calculate flux
        !--------------------------------------------------
        face%flux(1) = dt*sum(weight*vn*H_w*delta)+T1*sum(weight*vn*H_g*(1-delta))+T4*sum(weight*vn*h*(1-delta))-T5*sum(weight*vn**2*cell%sh(:,:,idx)*(1-delta))
        face%flux(2) = dt*sum(weight*vn*vn*H_w*delta)&
                        +T1*sum(weight*vn*vn*H_g*(1-delta))+T4*sum(weight*vn*vn*h*(1-delta))-T5*sum(weight*vn*vn**2*cell%sh(:,:,idx)*(1-delta))
        face%flux(3) = dt*sum(weight*vn*vt*H_w*delta)&
                        +T1*sum(weight*vn*vt*H_g*(1-delta))+T4*sum(weight*vt*vn*h*(1-delta))-T5*sum(weight*vt*vn**2*cell%sh(:,:,idx)*(1-delta))
        face%flux(4) = dt*0.5*sum(weight*vn*((vn**2+vt**2)*H_w+B_w)*delta)+T1*0.5*sum(weight*vn*((vn**2+vt**2)*H_g+B_g)*(1-delta))&
                        +T4*0.5*sum(weight*vn*((vn**2+vt**2)*h+b)*(1-delta))&
                        -T5*0.5*sum(weight*vn**2*((vn**2+vt**2)*cell%sh(:,:,idx)+cell%sb(:,:,idx))*(1-delta))

        face%flux_h = dt*vn*H_w*delta+T1*vn*H_g*(1-delta)+T4*vn*h*(1-delta)-T5*vn**2*cell%sh(:,:,idx)*(1-delta)
        face%flux_b = dt*vn*B_w*delta+T1*vn*B_g*(1-delta)+T4*vn*b*(1-delta)-T5*vn**2*cell%sb(:,:,idx)*(1-delta)

        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !Convert to global frame
        
        face%flux = GlobalFrame(face%flux,face%cosx,face%cosy)

        !Total flux
        face%flux = face%length*face%flux
        face%flux_h = face%length*face%flux_h
        face%flux_b = face%length*face%flux_b
        
        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Deallocate array
        deallocate(vn)
        deallocate(vt)
        deallocate(delta)
        deallocate(h)
        deallocate(b)
        deallocate(H_g)
        deallocate(B_g)
        deallocate(H_w)
        deallocate(B_w)
    end subroutine MultiscaleFluxBoundary

    !--------------------------------------------------
    !>Calculate micro slope of Maxwellian distribution
    !>@param[in] prim        :primary variables
    !>@param[in] sw          :slope of conVars
    !>@return    MicroSlope  :slope of Maxwellian distribution
    !--------------------------------------------------
    function MicroSlope(prim,sw)
        real(KREAL), intent(in)                         :: prim(4),sw(4)
        real(KREAL)                                     :: MicroSlope(4)

        MicroSlope(4) = 4.0*prim(4)**2/((CK+2)*prim(1))*(2.0*sw(4)-2.0*prim(2)*sw(2)-2.0*prim(3)*sw(3)+sw(1)*(prim(2)**2+prim(3)**2-0.5*(CK+2)/prim(4)))
        MicroSlope(3) = 2.0*prim(4)/prim(1)*(sw(3)-prim(3)*sw(1))-prim(3)*MicroSlope(4)
        MicroSlope(2) = 2.0*prim(4)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*MicroSlope(4)
        MicroSlope(1) = sw(1)/prim(1)-prim(2)*MicroSlope(2)-prim(3)*MicroSlope(3)-0.5*(prim(2)**2+prim(3)**2+0.5*(CK+2)/prim(4))*MicroSlope(4)
    end function MicroSlope

    !--------------------------------------------------
    !>Calculate collision time
    !>@param[in] prim    :primary variables
    !>@return    GetTau  :collision time
    !--------------------------------------------------
    function GetTau(prim)
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL)                                     :: GetTau

        GetTau = MU_REF*2.0*prim(4)**(1-OMEGA)/prim(1)
    end function GetTau
    
    !--------------------------------------------------
    !>calculate moments of velocity and \xi
    !>@param[in] prim       :primary variables
    !>@param[out] Mu,Mv     :<u^n>,<v^m>
    !>@param[out] Mxi       :<\xi^2n>
    !>@param[out] MuL,MuR   :<u^n>_{>0},<u^n>_{<0}
    !--------------------------------------------------
    subroutine CalcMoment(prim,Mu,Mv,Mxi,MuL,MuR)
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL), intent(out)                        :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM)
        real(KREAL), intent(out)                        :: Mv(0:MTUM)
        real(KREAL), intent(out)                        :: Mxi(0:2)
        integer :: i

        !Moments of normal velocity
        MuL(0) = 0.5*erfc(-sqrt(prim(4))*prim(2))
        MuL(1) = prim(2)*MuL(0)+0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))
        MuR(0) = 0.5*erfc(sqrt(prim(4))*prim(2))
        MuR(1) = prim(2)*MuR(0)-0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))

        do i=2,MNUM
            MuL(i) = prim(2)*MuL(i-1)+0.5*(i-1)*MuL(i-2)/prim(4)
            MuR(i) = prim(2)*MuR(i-1)+0.5*(i-1)*MuR(i-2)/prim(4)
        end do

        Mu = MuL+MuR

        !Moments of tangential velocity
        Mv(0) = 1.0
        Mv(1) = prim(3)

        do i=2,MTUM
            Mv(i) = prim(3)*Mv(i-1)+0.5*(i-1)*Mv(i-2)/prim(4)
        end do

        !Moments of \xi
        Mxi(0) = 1.0 !<\xi^0>
        Mxi(1) = 0.5*CK/prim(4) !<\xi^2>
        Mxi(2) = CK*(CK+2.0)/(4.0*prim(4)**2) !<\xi^4>
    end subroutine CalcMoment

    !--------------------------------------------------
    !>Calculate <u^\alpha*v^\beta*\xi^\delta*\psi>
    !>@param[in] Mu,Mv      :<u^\alpha>,<v^\beta>
    !>@param[in] Mxi        :<\xi^delta>
    !>@param[in] alpha,beta :exponential index of u and v
    !>@param[in] delta      :exponential index of \xi
    !>@return    Moment_uvxi :moment of <u^\alpha*v^\beta*\xi^\delta*\psi>
    !--------------------------------------------------
    function Moment_uvxi(Mu,Mv,Mxi,alpha,beta,delta)
        real(KREAL), intent(in)                         :: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
        integer(KINT), intent(in)                       :: alpha,beta,delta
        real(KREAL)                                     :: Moment_uvxi(4)

        Moment_uvxi(1) = Mu(alpha)*Mv(beta)*Mxi(delta/2)
        Moment_uvxi(2) = Mu(alpha+1)*Mv(beta)*Mxi(delta/2)
        Moment_uvxi(3) = Mu(alpha)*Mv(beta+1)*Mxi(delta/2)
        Moment_uvxi(4) = 0.5*(Mu(alpha+2)*Mv(beta)*Mxi(delta/2)+Mu(alpha)*Mv(beta+2)*Mxi(delta/2)+Mu(alpha)*Mv(beta)*Mxi((delta+2)/2))
    end function Moment_uvxi

    !--------------------------------------------------
    !>Calculate <a*u^\alpha*v^\beta>
    !>@param[in] a          :micro slope of Maxwellian
    !>@param[in] Mu,Mv      :<u^\alpha>,<v^\beta>
    !>@param[in] Mxi        :<\xi^l>
    !>@param[in] alpha,beta :exponential index of u and v
    !>@return    Moment_auvxi  :moment of <a*u^\alpha*v^\beta>
    !--------------------------------------------------
    function Moment_auvxi(a,Mu,Mv,Mxi,alpha,beta)
        real(KREAL), intent(in)                         :: a(4)
        real(KREAL), intent(in)                         :: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
        integer(KINT), intent(in)                       :: alpha,beta
        real(KREAL)                                     :: Moment_auvxi(4)

        Moment_auvxi = a(1)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+0,0)+&
                    a(2)*Moment_uvxi(Mu,Mv,Mxi,alpha+1,beta+0,0)+&
                    a(3)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+1,0)+&
                    0.5*a(4)*Moment_uvxi(Mu,Mv,Mxi,alpha+2,beta+0,0)+&
                    0.5*a(4)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+2,0)+&
                    0.5*a(4)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+0,2)
    end function Moment_auvxi
end module Flux

!--------------------------------------------------
!>UGKS solver
!--------------------------------------------------
module Solver
    use Flux
    implicit none
    real(KREAL)                                         :: old_sumRes(4) = 0.0

contains
    !--------------------------------------------------
    !>Interpolation of the boundary cell
    !>@param[inout] targetCell    :the target boundary cell
    !>@param[inout] leftCell      :the left cell
    !>@param[inout] rightCell     :the right cell
    !>@param[in]    idx           :the index indicating i or j direction
    !--------------------------------------------------
    subroutine InterpBoundary(leftCell,targetCell,rightCell,idx,rot,bc)
        type(CellCenter), intent(inout)                 :: targetCell
        type(CellCenter), intent(inout)                 :: leftCell,rightCell
        integer(KINT), intent(in)                       :: idx
        real(KREAL), intent(in)                         :: rot
        real(KREAL), intent(in)                         :: bc(4) !Primary variables at boundary
        real(KREAL)                                     :: prim(4),tau
        real(KREAL), allocatable, dimension(:,:)        :: H0,B0 !Maxwellian distribution function

        
        if (RECONSTRUCTION_METHOD==LIMITER) then
            targetCell%sh(:,:,idx) = (rightCell%h-leftCell%h)/(0.5*rightCell%length(idx)+0.5*leftCell%length(idx))
            targetCell%sb(:,:,idx) = (rightCell%b-leftCell%b)/(0.5*rightCell%length(idx)+0.5*leftCell%length(idx))
        elseif (RECONSTRUCTION_METHOD==CENTRAL) then
            targetCell%sh(:,:,idx) = (rightCell%h-leftCell%h)/(0.5*rightCell%length(idx)+0.5*leftCell%length(idx))
            targetCell%sb(:,:,idx) = (rightCell%b-leftCell%b)/(0.5*rightCell%length(idx)+0.5*leftCell%length(idx))
        else
            stop "Error in RECONSTRUCTION_METHOD!"
        end if
    end subroutine InterpBoundary

    !--------------------------------------------------
    !>Interpolation of the inner cells
    !>@param[in]    leftCell      :the left cell
    !>@param[inout] targetCell    :the target cell
    !>@param[in]    rightCell     :the right cell
    !>@param[in]    idx           :the index indicating i or j direction
    !--------------------------------------------------
    subroutine InterpInner(leftCell,targetCell,rightCell,idx)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellCenter), intent(inout)                 :: targetCell
        integer(KINT), intent(in)                       :: idx

        if (RECONSTRUCTION_METHOD==LIMITER) then
            call VanLeerLimiter(leftCell,targetCell,rightCell,idx)
        elseif (RECONSTRUCTION_METHOD==CENTRAL) then
            targetCell%sh(:,:,idx) = (rightCell%h-leftCell%h)/(0.5*rightCell%length(idx)+targetCell%length(idx)+0.5*leftCell%length(idx))
            targetCell%sb(:,:,idx) = (rightCell%b-leftCell%b)/(0.5*rightCell%length(idx)+targetCell%length(idx)+0.5*leftCell%length(idx))
        else
            stop "Error in RECONSTRUCTION_METHOD!"
        end if
    end subroutine InterpInner

    !--------------------------------------------------
    !>Reconstruct the slope of initial distribution function
    !>Index method
    !---------------------------------
    !           (i,j+1)              |
    !      ----------------          |
    !      |              |          |
    !      |              |          |
    !      |              |          |
    ! (i,j)|     (i,j)    |(i+1,j)   |
    !      |      Cell    |          |
    !      |              |          |
    !      |              |          |
    !      ----------------          |
    !            (i,j)               |
    !---------------------------------
    !--------------------------------------------------
    subroutine Reconstruction()
        integer(KINT)                                   :: i,j

        if (RECONSTRUCTION_METHOD==FIRST_ORDER) return

        !$omp parallel
        !--------------------------------------------------
        !i direction
        !--------------------------------------------------
        
        !Boundary part
        !$omp do
        do j=IYMIN,IYMAX
            call InterpBoundary(ctr(IXMIN,j),ctr(IXMIN,j),ctr(IXMIN+1,j),IDIRC,RN,BC_W)
            call InterpBoundary(ctr(IXMAX-1,j),ctr(IXMAX,j),ctr(IXMAX,j),IDIRC,RY,BC_E)
        end do
        !$omp end do nowait

        !Inner part
        !$omp do
        do j=IYMIN,IYMAX
            do i=IXMIN+1,IXMAX-1
                call InterpInner(ctr(i-1,j),ctr(i,j),ctr(i+1,j),IDIRC)
            end do
        end do
        !$omp end do nowait

        !--------------------------------------------------
        !j direction
        !--------------------------------------------------

        !Boundary part
        !$omp do
        do i=IXMIN,IXMAX
            call InterpBoundary(ctr(i,IYMIN),ctr(i,IYMIN),ctr(i,IYMIN+1),JDIRC,RN,BC_S)
            call InterpBoundary(ctr(i,IYMAX-1),ctr(i,IYMAX),ctr(i,IYMAX),JDIRC,RY,BC_N)
        end do
        !$omp end do nowait

        !$omp do
        do j=IYMIN+1,IYMAX-1
            do i=IXMIN,IXMAX
                call InterpInner(ctr(i,j-1),ctr(i,j),ctr(i,j+1),JDIRC)
            end do
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine Reconstruction

    !--------------------------------------------------
    !>Calculate the flux across the interfaces
    !--------------------------------------------------
    subroutine Evolution()
        integer(KINT)                                   :: i,j
        
        BC_N(2) = U0*cos(FRE*simTime) !North boundary
        !--------------------------------------------------
        !Calculate interface conVars for b_slope calculation
        !--------------------------------------------------
        if (FLUX_TYPE==MULTIDIMENSION) then
            !Inner part
            !$omp parallel
            !$omp do
            do j=IYMIN,IYMAX
                do i=IXMIN+1,IXMAX
                    call CalcFaceConvars(ctr(i-1,j),vface(i,j),ctr(i,j))
                end do
            end do
            !$omp end do nowait

            !Inner part
            !$omp do
            do j=IYMIN+1,IYMAX
                do i=IXMIN,IXMAX
                    call CalcFaceConvars(ctr(i,j-1),hface(i,j),ctr(i,j))
                end do
            end do
            !$omp end do nowait
            !$omp end parallel
        end if

        !--------------------------------------------------
        !Calculate interface flux
        !--------------------------------------------------
        if (FLUX_TYPE==MULTIDIMENSION) then
            !--------------------------------------------------
            !i direction
            !--------------------------------------------------
            !Boundary part
            !$omp parallel
            !$omp do
            do j=IYMIN,IYMAX
                call CalcFluxBoundary(BC_W,vface(IXMIN,j),ctr(IXMIN,j),IDIRC,RN) !RN means no frame rotation
                call CalcFluxBoundary(BC_E,vface(IXMAX+1,j),ctr(IXMAX,j),IDIRC,RY) !RY means with frame rotation
            end do
            !$omp end do nowait

            !Inner part
            !$omp do
            do j=IYMIN+1,IYMAX-1
                do i=IXMIN+1,IXMAX
                    call CalcFlux(ctr(i-1,j),vface(i,j),ctr(i,j),IDIRC,vface(i,j+1),vface(i,j-1),1) !idb=1, not boundary, full central differcence for b_slope
                end do
            end do
            !$omp end do nowait
            !$omp do
            do i=IXMIN+1,IXMAX
                call CalcFlux(ctr(i-1,IYMIN),vface(i,IYMIN),ctr(i,IYMIN),IDIRC,vface(i,IYMIN+1),vface(i,IYMIN),0) !idb=0, boundary, half central difference for b_slope
            end do
            !$omp end do nowait
            !$omp do
            do i=IXMIN+1,IXMAX
                call CalcFlux(ctr(i-1,IYMAX),vface(i,IYMAX),ctr(i,IYMAX),IDIRC,vface(i,IYMAX),vface(i,IYMAX-1),0) !idb=0, boundary, half central difference for b_slope
            end do
            !$omp end do nowait

            !--------------------------------------------------
            !j direction
            !--------------------------------------------------
            !Boundary part
            !$omp do
            do i=IXMIN,IXMAX
                call CalcFluxBoundary(BC_S,hface(i,IYMIN),ctr(i,IYMIN),JDIRC,RN) !RN means no frame rotation
                call CalcFluxBoundary(BC_N,hface(i,IYMAX+1),ctr(i,IYMAX),JDIRC,RY) !RY means with frame rotation 
            end do
            !$omp end do nowait

            !Inner part
            !$omp do
            do j=IYMIN+1,IYMAX
                do i=IXMIN+1,IXMAX-1
                    call CalcFlux(ctr(i,j-1),hface(i,j),ctr(i,j),JDIRC,hface(i-1,j),hface(i+1,j),1) !idb=1, not boundary, full central differcence for b_slope
                end do
            end do
            !$omp end do nowait
            !$omp do
            do j=IYMIN+1,IYMAX
                    call CalcFlux(ctr(IXMIN,j-1),hface(IXMIN,j),ctr(IXMIN,j),JDIRC,hface(IXMIN,j),hface(IXMIN+1,j),0) !idb=0, boundary, half central difference for b_slope
            end do
            !$omp end do nowait
            !$omp do
            do j=IYMIN+1,IYMAX
                    call CalcFlux(ctr(IXMAX,j-1),hface(IXMAX,j),ctr(IXMAX,j),JDIRC,hface(IXMAX-1,j),hface(IXMAX,j),0) !idb=0, boundary, half central difference for b_slope
            end do
            !$omp end do nowait
            !$omp end parallel
        elseif (FLUX_TYPE==SPLITTING) then
            !--------------------------------------------------
            !i direction
            !--------------------------------------------------
            
            !Boundary part
            !$omp parallel
            !$omp do
            do j=IYMIN,IYMAX
                call CalcFluxBoundary(BC_W,vface(IXMIN,j),ctr(IXMIN,j),IDIRC,RN) !RN means no frame rotation
                call CalcFluxBoundary(BC_E,vface(IXMAX+1,j),ctr(IXMAX,j),IDIRC,RY) !RY means with frame rotation
            end do
            !$omp end do nowait

            !Inner part
            !$omp do
            do j=IYMIN,IYMAX
                do i=IXMIN+1,IXMAX
                    call CalcFlux(ctr(i-1,j),vface(i,j),ctr(i,j),IDIRC)
                end do
            end do
            !$omp end do nowait

            !--------------------------------------------------
            !j direction
            !--------------------------------------------------

            !Boundary part
            !$omp do
            do i=IXMIN,IXMAX
                call CalcFluxBoundary(BC_S,hface(i,IYMIN),ctr(i,IYMIN),JDIRC,RN) !RN means no frame rotation
                call CalcFluxBoundary(BC_N,hface(i,IYMAX+1),ctr(i,IYMAX),JDIRC,RY) !RY means with frame rotation 
            end do
            !$omp end do nowait

            !Inner part
            !$omp do
            do j=IYMIN+1,IYMAX
                do i=IXMIN,IXMAX
                    call CalcFlux(ctr(i,j-1),hface(i,j),ctr(i,j),JDIRC)
                end do
            end do
            !$omp end do nowait
            !$omp end parallel
        else
            stop "Error in FLUX_TYPE!"
        end if
    end subroutine Evolution

    !--------------------------------------------------
    !>Update cell averaged values
    !--------------------------------------------------
    subroutine Update()
        real(KREAL), allocatable, dimension(:,:)        :: H_old,B_old !Equilibrium distribution at t=t^n
        real(KREAL), allocatable, dimension(:,:)        :: H,B !Equilibrium distribution at t=t^{n+1}
        real(KREAL), allocatable, dimension(:,:)        :: H_plus,B_plus !Shakhov part
        real(KREAL)                                     :: conVars_old(4),prim_old(4),prim(4) !Conversative and primary variables at t^n and t^{n+1}
        real(KREAL)                                     :: tau_old,tau !Collision time at t^n and t^{n+1}
        real(KREAL)                                     :: qf(2)
        real(KREAL)                                     :: sumRes(4),sumAvg(4)
        integer(KINT)                                   :: i,j,on

        !Allocate arrays
        allocate(H_old(uNum,vNum))
        allocate(B_old(uNum,vNum))
        allocate(H(uNum,vNum))
        allocate(B(uNum,vNum))
        allocate(H_plus(uNum,vNum))
        allocate(B_plus(uNum,vNum))

        !set initial value
        if (floor(simTime/TT)>num_period) then
            res = 0.0
            sumRes = 0.0
            sumAvg = 0.0
            on = 1
        else
            on = 0
        end if

        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX
                !--------------------------------------------------
                !Store conVars^n and calculate H^n,B^n,\tau^n
                !--------------------------------------------------
                conVars_old = ctr(i,j)%conVars !Store conVars^n
                prim_old = GetPrimary(conVars_old) !Convert to primary variables
                call DiscreteMaxwell(H_old,B_old,uSpace,vSpace,prim_old) !Calculate Maxwellian
                tau_old = GetTau(prim_old) !Calculate collision time \tau^n

                !--------------------------------------------------
                !Update conVars^{n+1} and Calculate H^{n+1},B^{n+1},\tau^{n+1}
                !--------------------------------------------------
                ctr(i,j)%conVars = ctr(i,j)%conVars+(vface(i,j)%flux-vface(i+1,j)%flux+hface(i,j)%flux-hface(i,j+1)%flux)/ctr(i,j)%area !Update conVars^{n+1}

                prim = GetPrimary(ctr(i,j)%conVars)
                call DiscreteMaxwell(H,B,uSpace,vSpace,prim)
                tau = GetTau(prim)

                !--------------------------------------------------
                !Record residual
                !--------------------------------------------------
                if (on==1) then
                    sumRes = sumRes+ctr(i,j)%conVars
                    sumAvg = sumAvg+abs(ctr(i,j)%conVars)
                end if

                !--------------------------------------------------
                !Shakhov part
                !--------------------------------------------------
                !Calculate heat flux at t=t^n
                qf = GetHeatFlux(ctr(i,j)%h,ctr(i,j)%b,uSpace,vSpace,prim_old) 

                !h^+ = H+H^+ at t=t^n
                call ShakhovPart(H_old,B_old,uSpace,vSpace,qf,prim_old,H_plus,B_plus) !H^+ and B^+
                H_old = H_old+H_plus !h^+
                B_old = B_old+B_plus !b^+

                !h^+ = H+H^+ at t=t^{n+1}
                call ShakhovPart(H,B,uSpace,vSpace,qf,prim,H_plus,B_plus)
                H = H+H_plus
                B = B+B_plus

                !--------------------------------------------------
                !Update distribution function
                !--------------------------------------------------
                ctr(i,j)%h = (ctr(i,j)%h+(vface(i,j)%flux_h-vface(i+1,j)%flux_h+hface(i,j)%flux_h-hface(i,j+1)%flux_h)/ctr(i,j)%area+&
                                dt*(H/tau*(1-exp(-dt/tau))+(H_old-ctr(i,j)%h)/tau_old*exp(-dt/tau)))/(1.0+dt/tau*(1-exp(-dt/tau)))
                ctr(i,j)%b = (ctr(i,j)%b+(vface(i,j)%flux_b-vface(i+1,j)%flux_b+hface(i,j)%flux_b-hface(i,j+1)%flux_b)/ctr(i,j)%area+&
                                dt*(B/tau*(1-exp(-dt/tau))+(B_old-ctr(i,j)%b)/tau_old*exp(-dt/tau)))/(1.0+dt/tau*(1-exp(-dt/tau)))
            end do
        end do

        !Calculate final residual
        if (on==1) then
            res = abs(sumRes-old_sumRes)/(sumAvg+SMV)
            old_sumRes = sumRes
            num_period = num_period+1
        end if

        !Deallocate arrays
        deallocate(H_old)
        deallocate(B_old)
        deallocate(H)
        deallocate(B)
        deallocate(H_plus)
        deallocate(B_plus)
    end subroutine Update
end module Solver

!--------------------------------------------------
!>Initialization of mesh and intial flow field
!--------------------------------------------------
module Initialization
    use ControlParameters
    use Tools
    implicit none

contains
    !--------------------------------------------------
    !>Main initialization subroutine
    !--------------------------------------------------
    subroutine Init()
        call InitMesh() !Initialize mesh
        call InitVelocity() !Initialize uSpace, vSpace and weights
        call InitAllocation(uNum,vNum) !Allocate discrete velocity space
        call InitFlowField() !Set the initial value
        dt = CFL*min(ctr(IXMIN,IYMIN)%length(1),ctr(IXMIN,IYMIN)%length(2))/max(U_MAX,V_MAX)
    end subroutine Init
    
    subroutine InitVelocity()
        if (QUADRATURE_TYPE==NEWTON) then
            call InitVelocityNewton(uNum,vNum) !Initialize discrete velocity space using Newton–Cotes formulas
        elseif (QUADRATURE_TYPE==GAUSS) then
            call InitVelocityGauss() !Initialize discrete velocity space using Gaussian-Hermite type quadrature
        elseif (QUADRATURE_TYPE==TRAPEZOID) then
            call InitVelocityTrapezoidal() !Initialize discrete velocity space using Trapezoidal rule
        else
            stop "Error in QUADRATURE_TYPE!"
        end if
    end subroutine InitVelocity

    !--------------------------------------------------
    !>Initialize mesh
    !--------------------------------------------------
    subroutine InitMesh()
        if (MESH_TYPE==UNIFORM) then
            call InitUniformMesh() !Initialize Uniform mesh
        elseif (MESH_TYPE==NONUNIFORM) then
            call InitNonUniformMesh()
        elseif (MESH_TYPE==RAND) then
            call InitRandomMesh()
        else
            stop "Error in MESH_TYPE!"
        end if
    end subroutine InitMesh

    !--------------------------------------------------
    !>Initialize uniform mesh
    !--------------------------------------------------
    subroutine InitUniformMesh()
        real(KREAL)                                     :: dx,dy
        integer(KINT)                                   :: i,j

        !Cell length
        dx = (X_END-X_START)/(IXMAX-IXMIN+1)
        dy = (Y_END-Y_START)/(IYMAX-IYMIN+1)

        !Geometry (node coordinate)
        forall(i=IXMIN:IXMAX+1,j=IYMIN:IYMAX+1)
                geometry(i,j)%x = X_START+(i-1)*dx
                geometry(i,j)%y = Y_START+(j-1)*dy
        end forall

        !Cell center
        forall(i=IXMIN:IXMAX,j=IYMIN:IYMAX)
            ctr(i,j)%x = X_START+(i-IXMIN+0.5)*dx
            ctr(i,j)%y = Y_START+(j-IYMIN+0.5)*dy
            ctr(i,j)%length(1) = dx
            ctr(i,j)%length(2) = dy
            ctr(i,j)%area = dx*dy
        end forall
        
        !Vertical interface
        forall(i=IXMIN:IXMAX+1,j=IYMIN:IYMAX)
            vface(i,j)%length = dy
            vface(i,j)%cosx = 1.0
            vface(i,j)%cosy = 0.0
        end forall

        !Horizontal interface
        forall(i=IXMIN:IXMAX,j=IYMIN:IYMAX+1)
            hface(i,j)%length = dx
            hface(i,j)%cosx = 0.0
            hface(i,j)%cosy = 1.0
        end forall
    end subroutine InitUniformMesh

    !--------------------------------------------------
    !>Initialize random mesh
    !--------------------------------------------------
    subroutine InitRandomMesh()
        real(KREAL)                                     :: dx(IXMIN:IXMAX),dy(IYMIN:IYMAX)
        integer(KINT)                                   :: i,j

        call random_number(dx)
        call random_number(dy)

        !Cell length
        dx = dx/sum(dx)*(X_END-X_START)
        dy = dy/sum(dy)*(Y_END-Y_START)

        !Geometry (node coordinate)
        geometry(IXMIN,IYMIN)%x = X_START
        geometry(IXMIN,IYMIN)%y = Y_START
        do i=IXMIN+1,IXMAX+1
            geometry(i,IYMIN)%x = geometry(i-1,IYMIN)%x+dx(i-1)
            geometry(i,IYMIN)%y = Y_START
        end do
        do j=IYMIN+1,IYMAX+1
            geometry(IXMIN,j)%x = X_START
            geometry(IXMIN,j)%y = geometry(IXMIN,j-1)%y+dy(j-1)
        end do
        do j=IYMIN+1,IYMAX+1
            do i=IXMIN+1,IXMAX+1
                geometry(i,j)%x = geometry(i-1,j)%x+dx(i-1)
                geometry(i,j)%y = geometry(i,j-1)%y+dy(j-1)
            end do
        end do

        !Cell center
        ctr(IXMIN,IYMIN)%x = X_START+0.5*dx(IXMIN)
        ctr(IXMIN,IYMIN)%y = Y_START+0.5*dy(IYMIN)
        ctr(IXMIN,IYMIN)%length(1) = dx(IXMIN)
        ctr(IXMIN,IYMIN)%length(2) = dy(IYMIN)
        ctr(IXMIN,IYMIN)%area = dx(IXMIN)*dy(IYMIN)
        do i=IXMIN+1,IXMAX
            ctr(i,IYMIN)%x = ctr(i-1,IYMIN)%x+dx(i)
            ctr(i,IYMIN)%y = Y_START+0.5*dy(IYMIN)
            ctr(i,IYMIN)%length(1) = dx(i)
            ctr(i,IYMIN)%length(2) = dy(IYMIN)
            ctr(i,IYMIN)%area = dx(i)*dy(IYMIN)
        end do
        do j=IYMIN+1,IYMAX
            ctr(IXMIN,j)%x = X_START+0.5*dx(IXMIN)
            ctr(IXMIN,j)%y = ctr(IXMIN,j-1)%y+dy(j)
            ctr(IXMIN,j)%length(1) = dx(IXMIN)
            ctr(IXMIN,j)%length(2) = dy(j)
            ctr(IXMIN,j)%area = dx(IXMIN)*dy(j)
        end do
        do j=IYMIN+1,IYMAX
            do i=IXMIN+1,IXMAX
                ctr(i,j)%x = ctr(i-1,j)%x+dx(i)
                ctr(i,j)%y = ctr(i,j-1)%y+dy(j)
                ctr(i,j)%length(1) = dx(i)
                ctr(i,j)%length(2) = dx(j)
                ctr(i,j)%area = dx(i)*dy(j)
            end do
        end do

        !Vertical interface
        forall(i=IXMIN:IXMAX+1,j=IYMIN:IYMAX)
            vface(i,j)%length = dy(j)
            vface(i,j)%cosx = 1.0
            vface(i,j)%cosy = 0.0
        end forall

        !Horizontal interface
        forall(i=IXMIN:IXMAX,j=IYMIN:IYMAX+1)
            hface(i,j)%length = dx(i)
            hface(i,j)%cosx = 0.0
            hface(i,j)%cosy = 1.0
        end forall
    end subroutine InitRandomMesh

    !--------------------------------------------------
    !>Initialize Nonuniform mesh
    !--------------------------------------------------
    subroutine InitNonUniformMesh()
        real(KREAL)                                     :: dx(IXMIN:IXMAX),dy(IYMIN:IYMAX)
        real(KREAL)                                     :: x(IXMIN:IXMAX+1),y(IYMIN:IYMAX+1)
        real(KREAL)                                     :: ax,ay !The larger the a, the smaller the mesh size near the walls.
        integer(KINT)                                   :: i,j

        ax = 2.0
        ay = 3.5

        !Cell length
        x = (/(i,i=IXMIN-1,IXMAX)/)
        y = (/(j,j=IYMIN-1,IYMAX)/)
        x = x/(IXMAX-IXMIN+1)
        y = y/(IYMAX-IYMIN+1)
        x = 0.5+0.5*tanh(ax*(x-0.5))/tanh(ax*0.5)
        y = 0.5+0.5*tanh(ay*(y-0.5))/tanh(ay*0.5)
        do i=IXMIN,IXMAX
            dx(i) = x(i+1)-x(i)
        end do
        do j=IYMIN,IYMAX
            dy(j) = y(j+1)-y(j)
        end do

        !Geometry (node coordinate)
        forall(i=IXMIN:IXMAX+1,j=IYMIN:IYMAX+1)
                geometry(i,j)%x = x(i)
                geometry(i,j)%y = y(j)
        end forall 

        !Cell center
        ctr(IXMIN,IYMIN)%x = X_START+0.5*dx(IXMIN)
        ctr(IXMIN,IYMIN)%y = Y_START+0.5*dy(IYMIN)
        ctr(IXMIN,IYMIN)%length(1) = dx(IXMIN)
        ctr(IXMIN,IYMIN)%length(2) = dy(IYMIN)
        ctr(IXMIN,IYMIN)%area = dx(IXMIN)*dy(IYMIN)
        do i=IXMIN+1,IXMAX
            ctr(i,IYMIN)%x = ctr(i-1,IYMIN)%x+dx(i)
            ctr(i,IYMIN)%y = Y_START+0.5*dy(IYMIN)
            ctr(i,IYMIN)%length(1) = dx(i)
            ctr(i,IYMIN)%length(2) = dy(IYMIN)
            ctr(i,IYMIN)%area = dx(i)*dy(IYMIN)
        end do
        do j=IYMIN+1,IYMAX
            ctr(IXMIN,j)%x = X_START+0.5*dx(IXMIN)
            ctr(IXMIN,j)%y = ctr(IXMIN,j-1)%y+dy(j)
            ctr(IXMIN,j)%length(1) = dx(IXMIN)
            ctr(IXMIN,j)%length(2) = dy(j)
            ctr(IXMIN,j)%area = dx(IXMIN)*dy(j)
        end do
        do j=IYMIN+1,IYMAX
            do i=IXMIN+1,IXMAX
                ctr(i,j)%x = ctr(i-1,j)%x+dx(i)
                ctr(i,j)%y = ctr(i,j-1)%y+dy(j)
                ctr(i,j)%length(1) = dx(i)
                ctr(i,j)%length(2) = dx(j)
                ctr(i,j)%area = dx(i)*dy(j)
            end do
        end do

        !Vertical interface
        forall(i=IXMIN:IXMAX+1,j=IYMIN:IYMAX)
            vface(i,j)%length = dy(j)
            vface(i,j)%cosx = 1.0
            vface(i,j)%cosy = 0.0
        end forall

        !Horizontal interface
        forall(i=IXMIN:IXMAX,j=IYMIN:IYMAX+1)
            hface(i,j)%length = dx(i)
            hface(i,j)%cosx = 0.0
            hface(i,j)%cosy = 1.0
        end forall
    end subroutine InitNonUniformMesh

    !--------------------------------------------------
    !>Initialize discrete velocity space using Newton–Cotes formulas
    !--------------------------------------------------
    subroutine InitVelocityNewton(num_u,num_v)
        integer(KINT), intent(inout)                    :: num_u,num_v
        real(KREAL)                                     :: du,dv !Spacing in u and v velocity space
        integer(KINT)                                   :: i,j

        !Modify num_u and num_v if not appropriate
        num_u = (num_u/4)*4+1
        num_v = (num_v/4)*4+1

        !Allocate array
        allocate(uSpace(num_u,num_v))
        allocate(vSpace(num_u,num_v))
        allocate(weight(num_u,num_v))

        !spacing in u and v velocity space
        du = (U_MAX-U_MIN)/(num_u-1)
        dv = (V_MAX-V_MIN)/(num_v-1)

        !velocity space
        forall(i=1:num_u,j=1:num_v)
            uSpace(i,j) = U_MIN+(i-1)*du
            vSpace(i,j) = V_MIN+(j-1)*dv
            weight(i,j) = (newtonCoeff(i,num_u)*du)*(newtonCoeff(j,num_v)*dv)
        end forall

        contains
            !--------------------------------------------------
            !>Calculate the coefficient for newton-cotes formula
            !>@param[in] idx          :index in velocity space
            !>@param[in] num          :total number in velocity space
            !>@return    newtonCoeff  :coefficient for newton-cotes formula
            !--------------------------------------------------
            pure function newtonCoeff(idx,num)
                integer(KINT), intent(in)               :: idx,num
                real(KREAL)                             :: newtonCoeff

                if (idx==1 .or. idx==num) then 
                    newtonCoeff = 14.0/45.0
                else if (mod(idx-5,4)==0) then
                    newtonCoeff = 28.0/45.0
                else if (mod(idx-3,4)==0) then
                    newtonCoeff = 24.0/45.0
                else
                    newtonCoeff = 64.0/45.0
                end if
            end function newtonCoeff
    end subroutine InitVelocityNewton

    !--------------------------------------------------
    !>Set discrete velocity space using Gaussian-Hermite type quadrature
    !>@param[in] umid,vmid :middle value of the velocity space, zero or macroscopic velocity
    !--------------------------------------------------
    subroutine InitVelocityGauss()
        real(KREAL)                                     :: umid,vmid
        real(KREAL)                                     :: vcoords(16), weights(16) !Velocity points and weight for 28 points (symmetry)
        integer(KINT)                                   :: i,j

        !Set 28*28 velocity points and weight
        ! vcoords = [ -0.5392407922630E+01, -0.4628038787602E+01, -0.3997895360339E+01, -0.3438309154336E+01,&
        !             -0.2926155234545E+01, -0.2450765117455E+01, -0.2007226518418E+01, -0.1594180474269E+01,&
        !             -0.1213086106429E+01, -0.8681075880846E+00, -0.5662379126244E+00, -0.3172834649517E+00,&
        !             -0.1331473976273E+00, -0.2574593750171E-01, +0.2574593750171E-01, +0.1331473976273E+00,&
        !             +0.3172834649517E+00, +0.5662379126244E+00, +0.8681075880846E+00, +0.1213086106429E+01,&
        !             +0.1594180474269E+01, +0.2007226518418E+01, +0.2450765117455E+01, +0.2926155234545E+01,&
        !             +0.3438309154336E+01, +0.3997895360339E+01, +0.4628038787602E+01, +0.5392407922630E+01 ]

        ! weights = [ +0.2070921821819E-12, +0.3391774320172E-09, +0.6744233894962E-07, +0.3916031412192E-05,&
        !             +0.9416408715712E-04, +0.1130613659204E-02, +0.7620883072174E-02, +0.3130804321888E-01,&
        !             +0.8355201801999E-01, +0.1528864568113E+00, +0.2012086859914E+00, +0.1976903952423E+00,&
        !             +0.1450007948865E+00, +0.6573088665062E-01, +0.6573088665062E-01, +0.1450007948865E+00,&
        !             +0.1976903952423E+00, +0.2012086859914E+00, +0.1528864568113E+00, +0.8355201801999E-01,&
        !             +0.3130804321888E-01, +0.7620883072174E-02, +0.1130613659204E-02, +0.9416408715712E-04,&
        !             +0.3916031412192E-05, +0.6744233894962E-07, +0.3391774320172E-09, +0.2070921821819E-12 ]

        !Set 16*16 velocity points and weight
        vcoords = [ -0.3686007162724397E+1, -0.2863133883708075E+1, -0.2183921153095858E+1, -0.1588855862270055E+1,&
                    -0.1064246312116224E+1, -0.6163028841823999, -0.2673983721677653, -0.5297864393185113E-1,&
                    0.5297864393185113E-1, 0.2673983721677653, 0.6163028841823999, 0.1064246312116224E+1,&
                    0.1588855862270055E+1, 0.2183921153095858E+1, 0.2863133883708075E+1, 0.3686007162724397E+1]

        weights = [ 0.1192596926595344E-5, 0.2020636491324107E-3, 0.5367935756025333E-2, 0.4481410991746290E-1,&
                    0.1574482826187903, 0.2759533979884218, 0.2683307544726388, 0.1341091884533595,&
                    0.1341091884533595, 0.2683307544726388, 0.2759533979884218, 0.1574482826187903,&
                    0.4481410991746290E-1, 0.5367935756025333E-2, 0.2020636491324107E-3, 0.1192596926595344E-5]

        !set grid number for u-velocity and v-velocity
        uNum = 16
        vNum = 16

        !allocate discrete velocity space
        allocate(uSpace(uNum,vNum)) !x direction
        allocate(vSpace(uNum,vNum)) !y direction
        allocate(weight(uNum,vNum)) !weight at u_k and v_l

        umid = 0.0
        vmid = 0.0
        !set velocity space and weight
        forall(i=1:uNum,j=1:vNum)
            uSpace(i,j) = umid+vcoords(i)
            vSpace(i,j) = vmid+vcoords(j)
            weight(i,j) = (weights(i)*exp(vcoords(i)**2))*(weights(j)*exp(vcoords(j)**2))
        end forall

        !store the maximum micro velocity
        U_MAX = maxval(abs(uSpace(:,1)))
        V_MAX = maxval(abs(vSpace(1,:)))
    end subroutine InitVelocityGauss

    !--------------------------------------------------
    !>Initialize discrete velocity space using Trapezoidal rule
    !--------------------------------------------------
    subroutine InitVelocityTrapezoidal()
        real(KREAL)                                     :: du,dv !Spacing in u and v velocity space
        integer(KINT)                                   :: i,j

        !Allocate array
        allocate(uSpace(uNum,vNum))
        allocate(vSpace(uNum,vNum))
        allocate(weight(uNum,vNum))

        !spacing in u and v velocity space
        du = (U_MAX-U_MIN)/(uNum-1)
        dv = (V_MAX-V_MIN)/(vNum-1)

        !velocity space
        forall(i=1:uNum,j=1:vNum)
            uSpace(i,j) = U_MIN+(i-1)*du
            vSpace(i,j) = V_MIN+(j-1)*dv
            weight(i,j) = (traCoeff(i,uNum)*du)*(traCoeff(j,vNum)*dv)
        end forall

        contains
            !--------------------------------------------------
            !>Calculate the coefficient for Trapezoidal rule
            !>@param[in] idx          :index in velocity space
            !>@param[in] num          :total number in velocity space
            !>@return    traCoeff     :coefficient for newton-cotes formula
            !--------------------------------------------------
            pure function traCoeff(idx,num)
                integer(KINT), intent(in)               :: idx,num
                real(KREAL)                             :: traCoeff

                if (idx==1 .or. idx==num) then 
                    traCoeff = 0.5
                else
                    traCoeff = 1.0
                end if
            end function traCoeff
    end subroutine InitVelocityTrapezoidal
    
    !--------------------------------------------------
    !>Allocate arrays in velocity space
    !--------------------------------------------------
    subroutine InitAllocation(num_u,num_v)
        integer(KINT), intent(in)                       :: num_u,num_v
        integer(KINT)                                   :: i,j

        !Cell center
        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX
                allocate(ctr(i,j)%h(num_u,num_v))
                allocate(ctr(i,j)%b(num_u,num_v))
                allocate(ctr(i,j)%sh(num_u,num_v,2))
                allocate(ctr(i,j)%sb(num_u,num_v,2))
                ctr(i,j)%sh = 0.0
                ctr(i,j)%sb = 0.0
            end do
        end do

        !Cell interface
        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX+1
                allocate(vface(i,j)%flux_h(num_u,num_v))
                allocate(vface(i,j)%flux_b(num_u,num_v))
            end do
        end do

        do j=IYMIN,IYMAX+1
            do i=IXMIN,IXMAX
                allocate(hface(i,j)%flux_h(num_u,num_v))
                allocate(hface(i,j)%flux_b(num_u,num_v))
            end do
        end do
    end subroutine InitAllocation

    !--------------------------------------------------
    !>Set the initial condition
    !--------------------------------------------------
    subroutine InitFlowField()
        real(KREAL), allocatable, dimension(:,:)        :: H,B !Reduced distribution functions
        real(KREAL)                                     :: conVars(4) !Conservative variables
        integer(KINT)                                   :: i,j

        !Allocation
        allocate(H(uNum,vNum))
        allocate(B(uNum,vNum))

        !Get conservative variables and Maxwellian distribution function
        conVars = GetConserved(INIT_GAS)
        call DiscreteMaxwell(H,B,uSpace,vSpace,INIT_GAS)

        !Initialize field
        forall(i=IXMIN:IXMAX,j=IYMIN:IYMAX)
            ctr(i,j)%conVars = conVars
            ctr(i,j)%h = H
            ctr(i,j)%b = B
            ctr(i,j)%sh = 0.0
            ctr(i,j)%sb = 0.0
        end forall

        !Deallocation
        deallocate(H)
        deallocate(B)
    end subroutine InitFlowField
    
    !--------------------------------------------------
    !>Deallocate arrays aftermath
    !--------------------------------------------------
    subroutine AfterDeallocation
        integer(KINT)                                   :: i,j

        !Deallocate array
        deallocate(uSpace)
        deallocate(vSpace)
        deallocate(weight)

        !Cell center
        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX
                deallocate(ctr(i,j)%h)
                deallocate(ctr(i,j)%b)
                deallocate(ctr(i,j)%sh)
                deallocate(ctr(i,j)%sb)
            end do
        end do

        !Cell interface
        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX+1
                deallocate(vface(i,j)%flux_h)
                deallocate(vface(i,j)%flux_b)
            end do
        end do
        do j=IYMIN,IYMAX+1
            do i=IXMIN,IXMAX
                deallocate(hface(i,j)%flux_h)
                deallocate(hface(i,j)%flux_b)
            end do
        end do
    end subroutine AfterDeallocation
end module Initialization

module Writer
    use Tools
    implicit none
    character(len=8)                                    :: date
    character(len=10)                                   :: time
    character(len=100)                                  :: fileName
contains
    !--------------------------------------------------
    !>Write result
    !--------------------------------------------------
    subroutine Output()
        real(KREAL)                                     :: prim(4)
        real(KREAL), dimension(:,:,:), allocatable      :: solution
        integer(KINT)                                   :: i,j,mid
        character(len=20)                               :: str

        !--------------------------------------------------
        !Prepare solutions
        !--------------------------------------------------
        allocate(solution(7,IXMIN:IXMAX,IYMIN:IYMAX))
        
        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX
                prim = GetPrimary(ctr(i,j)%conVars)
                solution(1:3,i,j) = prim(1:3) !Density,u,v
                solution(4,i,j) = 1/prim(4) !Temperature
                solution(5,i,j) = GetPressure(ctr(i,j)%h,ctr(i,j)%b,uSpace,vSpace,prim) !Pressure
                solution(6:7,i,j) = GetHeatFlux(ctr(i,j)%h,ctr(i,j)%b,uSpace,vSpace,prim) !Heat flux
            end do
        end do

        !--------------------------------------------------
        !Write to file
        !--------------------------------------------------
        !Open result file and write header
        !Using keyword arguments
        write(str , *) iter

        !Open result file
        open(unit=RSTFILE,file=RSTFILENAME//trim(fileName)//'_'//trim(adjustl(str))//'.dat',status="replace",action="write")
        
        !Write header
        write(RSTFILE,*) "VARIABLES = X, Y, Density, U, V, T, P, QX, QY"

        select case(OUTPUT_METHOD)
            case(CENTER)
                write(RSTFILE,*) "ZONE  I = ",IXMAX-IXMIN+2,", J = ",IYMAX-IYMIN+2,"DATAPACKING=BLOCK, VARLOCATION=([3-9]=CELLCENTERED)"

                !write geometry (node value)
                write(RSTFILE,"(6(ES23.16,2X))") geometry%x
                write(RSTFILE,"(6(ES23.16,2X))") geometry%y
            case(POINTS)
                write(RSTFILE,*) "ZONE  I = ",IXMAX-IXMIN+1,", J = ",IYMAX-IYMIN+1,"DATAPACKING=BLOCK"

                !write geometry (cell centered value)
                write(RSTFILE,"(6(ES23.16,2X))") ctr%x
                write(RSTFILE,"(6(ES23.16,2X))") ctr%y
        end select

        !Write solution (cell-centered)
        do i=1,7
            write(RSTFILE,"(6(ES23.16,2X))") solution(i,:,:)
        end do

        !close file
        close(RSTFILE)

        !Write central line
        mid = (IXMIN+IXMAX)/2
        open(unit=1,file='Y-U.plt',status="replace",action="write")
            write(1,*) 0.0,0.0
            do j=IYMIN,IYMAX
                write(1,*)  ctr(mid,j)%y,solution(2,mid,j)/0.15
            end do
            write(1,*) 1.0,1.0

        mid = (IYMIN+IYMAX)/2
        open(unit=1,file='X-V.plt',status="replace",action="write")
            write(1,*) 0.0,0.0
            do i=IXMIN,IXMAX
                write(1,*)  ctr(i,mid)%x,solution(3,i,mid)/0.15
            end do
            write(1,*) 1.0,0.0
        close(1)

        deallocate(solution)
    end subroutine Output
end module Writer

!--------------------------------------------------
!>Main program
!--------------------------------------------------
program Oscillatory_Cavity
    use Initialization
    use Solver
    use Writer
    implicit none
    real(KREAL)                                         :: start, finish
    integer(KINT)                                       :: old_num_period = -1
    
    !Initialization
    call Init()

    !Open history file and write header
    call date_and_time(DATE=date,TIME=time)
    fileName = '_'//date//'_'//time(1:6)
    open(unit=HSTFILE,file=HSTFILENAME//trim(fileName)//'.hst',status="replace",action="write") !Open history file
    write(HSTFILE,*) "VARIABLES = iter, simTime, dt" !write header

    !Open residual file and write header
    open(unit=RESFILE,file=RESFILENAME//trim(fileName)//'_residual.dat',status="replace",action="write") !Open residual file
    write(RESFILE,*) "VARIABLES = iter, resRho, resU, resV, resLambda" !write header
    close(RESFILE)

    !Star timer
    call cpu_time(start)

    !Iteration
    do while(.true.)
        call Reconstruction() !Calculate the slope of distribution function
        call Evolution() !Calculate flux across the interfaces
        call Update() !Update cell averaged value

        !Check stopping criterion
        if(all(res<EPS) .or. iter>=MAX_ITER) exit

        !Log the iteration situation every 10 iterations
        if (old_num_period<num_period) then
            write(*,"(A28,I15,2E15.7,I8)") "iter,simTime,dt,num_period:",iter,simTime,dt,num_period
            write(*,"(A18,4E15.7)") "res:",res
            write(HSTFILE,"(I15,2E15.7,I8)") iter,simTime,dt,num_period

            !Output the residual curve
            open(unit=RESFILE,file=RESFILENAME//trim(fileName)//'_residual.dat',status="old",action="write",position="append") !Open residual file
            write(RESFILE,"(I15,4E15.7)") iter,res,num_period
            close(RESFILE)
            old_num_period = num_period
        end if

        if (num_period*TT<=simTime .and. simTime<num_period*TT+dt) then
            call Output()
        elseif (num_period*TT+0.25*TT<=simTime .and. simTime<num_period*TT+0.25*TT+dt) then
            call Output()
        elseif (num_period*TT+0.5*TT<=simTime .and. simTime<num_period*TT+0.5*TT+dt) then
            call Output()
        elseif (num_period*TT+0.75*TT<=simTime .and. simTime<num_period*TT+0.75*TT+dt) then
            call Output()
        end if

        iter = iter+1
        simTime = simTime+dt
    end do

    !End timer
    call cpu_time(finish)
    print '("Run Time = ",f20.3," seconds.")', finish-start

    !Close history file
    close(HSTFILE)

    !Output solution
    call Output()

    !Aftermath
    call AfterDeallocation()
end program Oscillatory_Cavity