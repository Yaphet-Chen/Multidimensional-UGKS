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
    integer(KINT), parameter                            :: NONUNIFORM = 1 !Non-uniform mesh

    !Output
    integer(KINT), parameter                            :: CENTER = 0 !Output solution as cell centered value
    integer(KINT), parameter                            :: POINTS = 1 !Output solution as point value

    !Quadrature method
    integer(KINT), parameter                            :: NEWTON = 0 !Newton–Cotes
    integer(KINT), parameter                            :: GAUSS = 1 !Gauss-Hermite

    !Direction
    integer(KINT), parameter                            :: IDIRC = 1 !I direction
    integer(KINT), parameter                            :: JDIRC = 2 !J direction

    !Rotation
    integer(KINT), parameter                            :: RN = 1 !No frame rotation
    integer(KINT), parameter                            :: RY = -1 !With frame rotation
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
    integer(KINT), parameter                            :: OUTPUT_METHOD = POINTS
    real(KREAL), parameter                              :: CFL = 0.8 !CFL number
    integer(KINT), parameter                            :: MAX_ITER = 5E8 !Maximal iteration number
    real(KREAL), parameter                              :: EPS = 1.0E-5 !Convergence criteria
    real(KREAL)                                         :: simTime = 0.0 !Current simulation time
    integer(KINT)                                       :: iter = 1 !Number of iteration
    real(KREAL)                                         :: dt !Global time step
    real(KREAL)                                         :: res(4) !Residual
    
    !Output control
    character(len=13), parameter                        :: HSTFILENAME = "BoundaryLayer" !History file name
    character(len=13), parameter                        :: RSTFILENAME = "BoundaryLayer" !Result file name
    character(len=13), parameter                        :: RESFILENAME = "BoundaryLayer" !Residual file name
    integer(KINT), parameter                            :: HSTFILE = 20 !History file ID
    integer(KINT), parameter                            :: RSTFILE = 21 !Result file ID
    integer(KINT), parameter                            :: RESFILE = 22 !Residual file ID

    !Gas propeties
    integer(KINT), parameter                            :: CK = 0 !Internal degree of freedom, here 1 denotes monatomic gas
    real(KREAL), parameter                              :: GAMMA = real(CK+4,KREAL)/real(CK+2,KREAL) !Ratio of specific heat
    real(KREAL), parameter                              :: OMEGA = 0.0 !Temperature dependence index in HS/VHS/VSS model
    real(KREAL), parameter                              :: PR = 2.0/3.0 !Prandtl number

    !MU_REF determined by Re number
    real(KREAL), parameter                              :: MA = 0.3 !Free stream inflow Mach number
    real(KREAL), parameter                              :: Re = 100000.0 !Reynolds number in reference state
    real(KREAL), parameter                              :: MU_REF = MA*sqrt(0.5*GAMMA)*104.68/Re !Viscosity coefficient in reference state

    !Geometry
    real(KREAL), parameter                              :: X_START = 0.0, Y_START = 0.0 !Start point in x, y direction
    real(KREAL), parameter                              :: RX_L = 1.1, RX_R = 1.05, RY_U = 1.2 !Common ratio at x and y direction
    real(KREAL), parameter                              :: DX_MIN = 0.1, DY_MIN = 0.0175 !Scale factor, i.e. minimal cell size
    integer(KINT), parameter                            :: X_NUM_L = 38, X_NUM_R = 82, Y_NUM = 30 !Points number in x, y direction
    integer(KINT), parameter                            :: IXMIN = -X_NUM_L+1 , IXMAX = X_NUM_R, IYMIN = 1 , IYMAX = Y_NUM !Cell index range
    integer(KINT), parameter                            :: GHOST = 2 !Ghost point number
    integer(KINT), parameter                            :: N_GRID = (IXMAX-IXMIN+1)*(IYMAX-IYMIN+1) !Total number of cell
    
    !--------------------------------------------------
    !Discrete velocity space
    !--------------------------------------------------
    integer(KINT)                                       :: uNum = 16, vNum = 16 !Number of points in velocity space for u and v
    real(KREAL)                                         :: U_MIN = -2.5, U_MAX = +2.5, V_MIN = -2.5, V_MAX = +2.5 !Minimum and maximum micro velocity
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
    type(CellCenter)                                    :: ctr(IXMIN-GHOST:IXMAX+GHOST,IYMIN-GHOST:IYMAX+GHOST) !Cell center
    type(CellInterface)                                 :: vface(IXMIN:IXMAX+1,IYMIN:IYMAX),hface(IXMIN:IXMAX,IYMIN:IYMAX+1) !Vertical and horizontal interfaces
    type(Grid)                                          :: geometry(IXMIN:IXMAX+1,IYMIN:IYMAX+1)
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
    integer(KREAL), parameter                           :: MNUM = 6 !Number of normal velocity moments
    integer(KREAL), parameter                           :: MTUM = 5 !Number of tangential velocity moments

contains
    !--------------------------------------------------
    !>Calculate flux of inner interface
    !>@param[in]    leftCell  :cell left to the target interface
    !>@param[inout] face      :the target interface
    !>@param[in]    rightCell :cell right to the target interface
    !>@param[in]    idx       :index indicating i or j direction
    !--------------------------------------------------
    subroutine CalcFlux(leftCell,face,rightCell,idx)
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
        real(KREAL)                                     :: a_slope(4),aT(4) !Micro slope of Maxwellian distribution, normal and time.
        real(KREAL)                                     :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM),Mv(0:MTUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<v^m>,<\xi^l>
        real(KREAL)                                     :: Mau0(4),Mau(4),MauT(4) !<u\psi>,<a*u^n*\psi>,<A*u*\psi>
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
        conVars = 0.5*LocalFrame(rightCell%conVars+leftCell%conVars,face%cosx,face%cosy)
        !Convert to primary variables
        prim = GetPrimary(conVars)

        !--------------------------------------------------
        !Calculate a_slope
        !--------------------------------------------------
        sw = LocalFrame(rightCell%conVars-leftCell%conVars,face%cosx,face%cosy)/(0.5*leftCell%length(idx)+0.5*rightCell%length(idx)) !left slope of conVars
        a_slope = MicroSlope(prim,sw) !Calculate a_slope

        !--------------------------------------------------
        !Calculate time slope of conVars and A
        !--------------------------------------------------
        !<u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
        call CalcMoment(prim,Mu,Mv,Mxi,MuL,MuR) 

        Mau = Moment_auvxi(a_slope,Mu,Mv,Mxi,1,0) !<a*u*\psi>

        sw = -prim(1)*Mau !Time slope of conVars
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
        Mau = Moment_auvxi(a_slope,Mu,Mv,Mxi,2,0) !<a*u^2*\psi>
        MauT = Moment_auvxi(aT,Mu,Mv,Mxi,1,0) !<A*u*\psi>

        face%flux = Mt(1)*prim(1)*Mau0+Mt(2)*prim(1)*Mau+Mt(3)*prim(1)*MauT

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
                        Mt(1)*0.5*sum(weight*vn*((vn**2+vt**2)*H_plus+B_plus))+&
                        Mt(4)*0.5*sum(weight*vn*((vn**2+vt**2)*h+b))-&
                        Mt(5)*0.5*sum(weight*vn**2*((vn**2+vt**2)*sh+sb))

        !--------------------------------------------------
        !Calculate flux of distribution function
        !--------------------------------------------------
        face%flux_h = Mt(1)*vn*(H0+H_plus)+&
                        Mt(2)*vn**2*(a_slope(1)*H0+a_slope(2)*vn*H0+a_slope(3)*vt*H0+0.5*a_slope(4)*((vn**2+vt**2)*H0+B0))+&
                        Mt(3)*vn*(aT(1)*H0+aT(2)*vn*H0+aT(3)*vt*H0+0.5*aT(4)*((vn**2+vt**2)*H0+B0))+&
                        Mt(4)*vn*h-Mt(5)*vn**2*sh

        face%flux_b = Mt(1)*vn*(B0+B_plus)+&
                        Mt(2)*vn**2*(a_slope(1)*B0+a_slope(2)*vn*B0+a_slope(3)*vt*B0+0.5*a_slope(4)*((vn**2+vt**2)*B0+Mxi(2)*H0))+&
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
    end subroutine CalcFlux

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
    !>Calculate <a*u^\alpha*v^\beta*\psi>
    !>@param[in] a          :micro slope of Maxwellian
    !>@param[in] Mu,Mv      :<u^\alpha>,<v^\beta>
    !>@param[in] Mxi        :<\xi^l>
    !>@param[in] alpha,beta :exponential index of u and v
    !>@return    Moment_auvxi  :moment of <a*u^\alpha*v^\beta*\psi>
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

contains
    !--------------------------------------------------
    !>Calculate time step
    !--------------------------------------------------
    subroutine TimeStep()
        real(KREAL)                                     :: prim(4) !Primary variables
        real(KREAL)                                     :: sos !Speed of sound
        real(KREAL)                                     :: tMax !Max 1/dt allowed
        integer(KINT)                                   :: i,j
        
        !Set initial value
        tMax = 0.0

        !$omp parallel 
        !$omp do private(i,j,sos,prim) reduction(max:tmax)
        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX
                !Convert conservative variables to primary variables
                prim = GetPrimary(ctr(i,j)%conVars)

                !Get sound speed
                sos = GetSoundSpeed(prim)

                !Maximum velocity
                prim(2) = max(U_MAX,abs(prim(2)))+sos
                prim(3) = max(V_MAX,abs(prim(3)))+sos

                !Maximum 1/dt allowed
                tMax = max(tMax,(ctr(i,j)%length(2)*prim(2)+ctr(i,j)%length(1)*prim(3))/ctr(i,j)%area)
            end do
        end do
        !$omp end do
        !$omp end parallel
        
        !Time step
        dt = CFL/tMax
    end subroutine TimeStep

    subroutine Boundary()
        call OutBoundary() !Set free stream outflow boundary at right and up
        call BottomBoundary() !Set adiabatic slip/non-slip boundary at bottom
    end subroutine Boundary

    !--------------------------------------------------
    !>Set adiabatic slip(x<0) or non-slip(x>o) boundary at bottom
    !--------------------------------------------------
    subroutine BottomBoundary()
        integer(KINT)                                   :: i,j,l,m

        !$omp parallel	
        !$omp do
        do j=1,GHOST
            do i=1,IXMAX
                ctr(i,IYMIN-j)%conVars(1) = ctr(i,IYMIN+j-1)%conVars(1)
                ctr(i,IYMIN-j)%conVars(2) = -ctr(i,IYMIN+j-1)%conVars(2) !Reverse u velocity
                ctr(i,IYMIN-j)%conVars(3) = -ctr(i,IYMIN+j-1)%conVars(3) !Reverse v velocity
                ctr(i,IYMIN-j)%conVars(4) = ctr(i,IYMIN+j-1)%conVars(4)
                do m=1,vNum
                    do l=1,uNum
                        !Make the distribution function central symmetry
                        ctr(i,IYMIN-j)%h(l,m) = ctr(i,IYMIN+j-1)%h(uNum-l+1,vNum-m+1)
                        ctr(i,IYMIN-j)%b(l,m) = ctr(i,IYMIN+j-1)%b(uNum-l+1,vNum-m+1)
                        !Do the same thing for slope
                        ctr(i,IYMIN-j)%sh = 0.0
                        ctr(i,IYMIN-j)%sb = 0.0
                    end do
                end do
            end do
        end do
        !$omp end do nowait	
        	
        !$omp do
        do j=1,GHOST
            do i=IXMIN,0
                ctr(i,IYMIN-j)%conVars(1) = ctr(i,IYMIN+j-1)%conVars(1)
                ctr(i,IYMIN-j)%conVars(2) = ctr(i,IYMIN+j-1)%conVars(2)
                ctr(i,IYMIN-j)%conVars(3) = -ctr(i,IYMIN+j-1)%conVars(3) !Only reverse v velocity
                ctr(i,IYMIN-j)%conVars(4) = ctr(i,IYMIN+j-1)%conVars(4)
                do m=1,vNum
                    do l=1,uNum
                        !Make the distribution function u axial symmetry
                        ctr(i,IYMIN-j)%h(l,m) = ctr(i,IYMIN+j-1)%h(l,vNum-m+1)
                        ctr(i,IYMIN-j)%b(l,m) = ctr(i,IYMIN+j-1)%b(l,vNum-m+1)
                        !Do the same thing for slope
                        ctr(i,IYMIN-j)%sh = 0.0
                        ctr(i,IYMIN-j)%sb = 0.0
                    end do
                end do
            end do
        end do
        !$omp end do nowait	
        !$omp end parallel
    end subroutine BottomBoundary
    
    subroutine OutBoundary()
        integer(KINT)                                   :: i,j

        !$omp parallel	
        !$omp do
        !Set upper free stream outflow boundary
        do j=1,GHOST
            do i=IXMIN,IXMAX
                ctr(i,IYMAX+j)%conVars = ctr(i,IYMAX)%conVars
                ctr(i,IYMAX+j)%h = ctr(i,IYMAX)%h
                ctr(i,IYMAX+j)%b = ctr(i,IYMAX)%b
                ctr(i,IYMAX+j)%sh = ctr(i,IYMAX)%sh
                ctr(i,IYMAX+j)%sb = ctr(i,IYMAX)%sb
            end do
        end do
        !$omp end do nowait	
        	
        !$omp do
        !Set right free stream outflow boundary
        do i=1,GHOST
            do j=IYMIN,IYMAX
                ctr(IXMAX+i,j)%conVars = ctr(IXMAX,j)%conVars
                ctr(IXMAX+i,j)%h = ctr(IXMAX,j)%h
                ctr(IXMAX+i,j)%b = ctr(IXMAX,j)%b
                ctr(IXMAX+i,j)%sh = ctr(IXMAX,j)%sh
                ctr(IXMAX+i,j)%sb = ctr(IXMAX,j)%sb
            end do
        end do
        !$omp end do nowait	
        !$omp end parallel
    end subroutine OutBoundary

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
        !Inner part
        !$omp do
        do j=IYMIN,IYMAX
            do i=IXMIN-GHOST+1,IXMAX+GHOST-1
                call InterpInner(ctr(i-1,j),ctr(i,j),ctr(i+1,j),IDIRC)
            end do
        end do
        !$omp end do nowait

        !--------------------------------------------------
        !j direction
        !--------------------------------------------------
        !$omp do
        do j=IYMIN-GHOST+1,IYMAX+GHOST-1
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
        !--------------------------------------------------
        !Calculate interface flux
        !--------------------------------------------------
        !--------------------------------------------------
        !i direction
        !--------------------------------------------------
        !Inner part
        !$omp parallel
        !$omp do
        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX+1
                call CalcFlux(ctr(i-1,j),vface(i,j),ctr(i,j),IDIRC)
            end do
        end do
        !$omp end do nowait

        !--------------------------------------------------
        !j direction
        !--------------------------------------------------
        !Inner part
        !$omp do
        do j=IYMIN,IYMAX+1
            do i=IXMIN,IXMAX
                call CalcFlux(ctr(i,j-1),hface(i,j),ctr(i,j),JDIRC)
            end do
        end do
        !$omp end do nowait
        !$omp end parallel
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
        integer(KINT)                                   :: i,j

        !Allocate arrays
        allocate(H_old(uNum,vNum))
        allocate(B_old(uNum,vNum))
        allocate(H(uNum,vNum))
        allocate(B(uNum,vNum))
        allocate(H_plus(uNum,vNum))
        allocate(B_plus(uNum,vNum))

        !set initial value
        res = 0.0
        sumRes = 0.0
        sumAvg = 0.0

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
                sumRes = sumRes+((conVars_old-ctr(i,j)%conVars)*ctr(i,j)%area)**2
                sumAvg = sumAvg+abs(ctr(i,j)%conVars)*ctr(i,j)%area
            
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
                ! ctr(i,j)%h = (ctr(i,j)%h+(vface(i,j)%flux_h-vface(i+1,j)%flux_h+hface(i,j)%flux_h-hface(i,j+1)%flux_h)/ctr(i,j)%area+&
                !                     0.5*dt*(H/tau+(H_old-ctr(i,j)%h)/tau_old))/(1.0+0.5*dt/tau)
                ! ctr(i,j)%b = (ctr(i,j)%b+(vface(i,j)%flux_b-vface(i+1,j)%flux_b+hface(i,j)%flux_b-hface(i,j+1)%flux_b)/ctr(i,j)%area+&
                !                     0.5*dt*(B/tau+(B_old-ctr(i,j)%b)/tau_old))/(1.0+0.5*dt/tau)

                ctr(i,j)%h = (ctr(i,j)%h+(vface(i,j)%flux_h-vface(i+1,j)%flux_h+hface(i,j)%flux_h-hface(i,j+1)%flux_h)/ctr(i,j)%area+&
                                    dt*(H/tau*(1-exp(-dt/tau))+(H_old-ctr(i,j)%h)/tau_old*exp(-dt/tau)))/(1.0+dt/tau*(1-exp(-dt/tau)))
                ctr(i,j)%b = (ctr(i,j)%b+(vface(i,j)%flux_b-vface(i+1,j)%flux_b+hface(i,j)%flux_b-hface(i,j+1)%flux_b)/ctr(i,j)%area+&
                                    dt*(B/tau*(1-exp(-dt/tau))+(B_old-ctr(i,j)%b)/tau_old*exp(-dt/tau)))/(1.0+dt/tau*(1-exp(-dt/tau)))
            end do
        end do

        !Calculate final residual
        res = sqrt(N_GRID*sumRes)/(sumAvg+SMV)/dt
        
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
    end subroutine Init
    
    subroutine InitVelocity()
        if (QUADRATURE_TYPE==NEWTON) then
            call InitVelocityNewton(uNum,vNum) !Initialize discrete velocity space using Newton–Cotes formulas
        elseif (QUADRATURE_TYPE==GAUSS) then
            call InitVelocityGauss() !Initialize discrete velocity space using Gaussian-Hermite type quadrature
        else
            stop "Error in QUADRATURE_TYPE!"
        end if
    end subroutine InitVelocity

    !--------------------------------------------------
    !>Initialize mesh
    !--------------------------------------------------
    subroutine InitMesh()
        if (MESH_TYPE==NONUNIFORM) then
            call InitNonUniformMesh()
        else
            stop "Error in MESH_TYPE!"
        end if
    end subroutine InitMesh

    !--------------------------------------------------
    !>Initialize Nonuniform mesh
    !--------------------------------------------------
    subroutine InitNonUniformMesh()
        integer(KINT)                                   :: i,j

        !Cell center
        !Y direction
        ctr(:,IYMIN-GHOST)%length(2) = DY_MIN
        ctr(:,IYMIN-GHOST)%y = Y_START-(RY_U*ctr(:,IYMIN-GHOST)%Length(2)+0.5*ctr(:,IYMIN-GHOST)%Length(2))
        do j=IYMIN-GHOST+1,IYMAX+GHOST
            ctr(:,j)%length(2) = ctr(:,j-1)%length(2)*RY_U
            ctr(:,j)%y = ctr(:,j-1)%y+0.5*(ctr(:,j)%length(2)+ctr(:,j-1)%length(2))
        end do
        
        !X direction (x>0)
        ctr(1,:)%x = X_START+0.5*DX_MIN
        ctr(1,:)%length(1) = DX_MIN
        do i=2,IXMAX+GHOST
            ctr(i,:)%length(1) = ctr(i-1,:)%length(1)*RX_R
            ctr(i,:)%x = ctr(i-1,:)%x+0.5*(ctr(i,:)%length(1)+ctr(i-1,:)%length(1))
        end do

        !X direction (x<0)
        do i=0,IXMIN-GHOST,-1
            ctr(i,:)%length(1) = ctr(i+1,:)%length(1)*RX_L
            ctr(i,:)%x = ctr(i+1,:)%x-0.5*(ctr(i,:)%length(1)+ctr(i+1,:)%length(1))
        end do

        do j=IYMIN-GHOST,IYMAX+GHOST
            do i=IXMIN-GHOST,IXMAX+GHOST
                ctr(i,j)%area = ctr(i,j)%length(1)*ctr(i,j)%length(2)
            end do
        end do

        !Vertical interface
        forall(i=IXMIN:IXMAX+1,j=IYMIN:IYMAX)
            vface(i,j)%length = ctr(i,j)%length(2)
            vface(i,j)%cosx = 1.0
            vface(i,j)%cosy = 0.0
        end forall

        !Horizontal interface
        forall(i=IXMIN:IXMAX,j=IYMIN:IYMAX+1)
            hface(i,j)%length = ctr(i,j)%length(1)
            hface(i,j)%cosx = 0.0
            hface(i,j)%cosy = 1.0
        end forall

        !Geometry (node coordinate)
        forall(i=IXMIN:IXMAX+1,j=IYMIN:IYMAX+1)
            geometry(i,j)%y = 0.5*(ctr(i,j-1)%y+ctr(i,j)%y)
            geometry(i,j)%x = 0.5*(ctr(i-1,j)%x+ctr(i,j)%x)
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

        !Set 28x28 velocity points and weight
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

        !Set 16x16 velocity points and weight
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
    !>Allocate arrays in velocity space
    !--------------------------------------------------
    subroutine InitAllocation(num_u,num_v)
        integer(KINT), intent(in)                       :: num_u,num_v
        integer(KINT)                                   :: i,j

        !Cell center
        do j=IYMIN-GHOST,IYMAX+GHOST
            do i=IXMIN-GHOST,IXMAX+GHOST
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
        real(KREAL), dimension(4)                       :: INIT_GAS = [1.0, 0.0, 0.0, 1.0] !Initial condition (density, u-velocity, v-velocity, lambda=1/temperature)
        integer(KINT)                                   :: i,j

        !Allocation
        allocate(H(uNum,vNum))
        allocate(B(uNum,vNum))

        !Get conservative variables and Maxwellian distribution function
        INIT_GAS(2) = MA*sqrt(0.5*GAMMA) !Set u-velocity
        conVars = GetConserved(INIT_GAS)
        call DiscreteMaxwell(H,B,uSpace,vSpace,INIT_GAS)

        !Initialize field and left free stream inflow boundary
        forall(i=IXMIN-GHOST:IXMAX+GHOST,j=IYMIN-GHOST:IYMAX+GHOST)
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
        do j=IYMIN-GHOST,IYMAX+GHOST
            do i=IXMIN-GHOST,IXMAX+GHOST
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
        integer(KINT)                                   :: i,j
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
                solution(5,i,j) = 0.5*solution(4,i,j)*solution(1,i,j) !Pressure
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
                write(RSTFILE,"(6(ES23.16,2X))") ctr(IXMIN:IXMAX,IYMIN:IYMAX)%x
                write(RSTFILE,"(6(ES23.16,2X))") ctr(IXMIN:IXMAX,IYMIN:IYMAX)%y
        end select

        !Write solution (cell-centered)
        do i=1,7
            write(RSTFILE,"(6(ES23.16,2X))") solution(i,:,:)
        end do

        !close file
        close(RSTFILE)

        open(unit=1,file='10.plt',status="replace",action="write")
            do j=IYMIN,IYMAX
                do i=IXMIN,IXMAX
                    if(i==10)then
                        write(1,*)  ctr(i,j)%y/sqrt(MU_REF*ctr(i,j)%x/(MA*sqrt(0.5*GAMMA))),solution(2,i,j)/(MA*sqrt(0.5*GAMMA)),solution(3,i,j)*sqrt(ctr(i,j)%x/MU_REF/(MA*sqrt(0.5*GAMMA)))
                    end if
                end do
            end do
        close(1)
        open(unit=1,file='20.plt',status="replace",action="write")
            do j=IYMIN,IYMAX
                do i=IXMIN,IXMAX
                    if(i==20)then
                        write(1,*)  ctr(i,j)%y/sqrt(MU_REF*ctr(i,j)%x/(MA*sqrt(0.5*GAMMA))),solution(2,i,j)/(MA*sqrt(0.5*GAMMA)),solution(3,i,j)*sqrt(ctr(i,j)%x/MU_REF/(MA*sqrt(0.5*GAMMA)))
                    end if
                end do
            end do
        close(1)
        open(unit=1,file='30.plt',status="replace",action="write")
            do j=IYMIN,IYMAX
                do i=IXMIN,IXMAX
                    if(i==30)then
                        write(1,*)  ctr(i,j)%y/sqrt(MU_REF*ctr(i,j)%x/(MA*sqrt(0.5*GAMMA))),solution(2,i,j)/(MA*sqrt(0.5*GAMMA)),solution(3,i,j)*sqrt(ctr(i,j)%x/MU_REF/(MA*sqrt(0.5*GAMMA)))
                    end if
                end do
            end do
        close(1)
        open(unit=1,file='40.plt',status="replace",action="write")
            do j=IYMIN,IYMAX
                do i=IXMIN,IXMAX
                    if(i==40)then
                        write(1,*)  ctr(i,j)%y/sqrt(MU_REF*ctr(i,j)%x/(MA*sqrt(0.5*GAMMA))),solution(2,i,j)/(MA*sqrt(0.5*GAMMA)),solution(3,i,j)*sqrt(ctr(i,j)%x/MU_REF/(MA*sqrt(0.5*GAMMA)))
                    end if
                end do
            end do
        close(1)
        open(unit=1,file='50.plt',status="replace",action="write")
            do j=IYMIN,IYMAX
                do i=IXMIN,IXMAX
                    if(i==50)then
                        write(1,*)  ctr(i,j)%y/sqrt(MU_REF*ctr(i,j)%x/(MA*sqrt(0.5*GAMMA))),solution(2,i,j)/(MA*sqrt(0.5*GAMMA)),solution(3,i,j)*sqrt(ctr(i,j)%x/MU_REF/(MA*sqrt(0.5*GAMMA)))
                    end if
                end do
            end do
        close(1)
        open(unit=1,file='60.plt',status="replace",action="write")
            do j=IYMIN,IYMAX
                do i=IXMIN,IXMAX
                    if(i==60)then
                        write(1,*)  ctr(i,j)%y/sqrt(MU_REF*ctr(i,j)%x/(MA*sqrt(0.5*GAMMA))),solution(2,i,j)/(MA*sqrt(0.5*GAMMA)),solution(3,i,j)*sqrt(ctr(i,j)%x/MU_REF/(MA*sqrt(0.5*GAMMA)))
                    end if
                end do
            end do
        close(1)

        deallocate(solution)
    end subroutine Output
end module Writer

!--------------------------------------------------
!>Main program
!--------------------------------------------------
program BoundaryLayer
    use Initialization
    use Solver
    use Writer
    implicit none
    real(KREAL)                                         :: start, finish
    
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
        call TimeStep() !Calculate the time step
        call Boundary() !Set Boundary condition
        call Reconstruction() !Calculate the slope of distribution function
        call Evolution() !Calculate flux across the interfaces
        call Update() !Update cell averaged value

        !Check stopping criterion
        if(all(res<EPS) .or. iter>=MAX_ITER) exit

        !Log the iteration situation every 10 iterations
        if (mod(iter,10)==0) then
            write(*,"(A18,I15,2E15.7)") "iter,simTime,dt:",iter,simTime,dt
            write(*,"(A18,4E15.7)") "res:",res
            write(HSTFILE,"(I15,2E15.7)") iter,simTime,dt

            !Output the residual curve
            open(unit=RESFILE,file=RESFILENAME//trim(fileName)//'_residual.dat',status="old",action="write",position="append") !Open residual file
            write(RESFILE,"(I15,4E15.7)") iter,res
            close(RESFILE)
        end if

        if (mod(iter,10000)==0) then
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
end program BoundaryLayer
