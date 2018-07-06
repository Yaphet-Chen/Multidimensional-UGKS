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

    !Direction
    integer(KINT), parameter                            :: IDIRC = 1 !I direction
    integer(KINT), parameter                            :: JDIRC = 2 !J direction
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
        !Flux
        real(KREAL)                                     :: flux(4) !Conservative variables flux at cell interface: density flux, x and y momentum flux, total energy flux
        real(KREAL), allocatable, dimension(:,:)        :: flux_h,flux_b !Flux of distribution function
    end type CellInterface

    !index method
    !          (i,j+1)
    !     ----------------
    !     |              |
    !     |              |
    !     |              |
    !(i,j)|     (i,j)    |(i+1,j)
    !     |      Cell    |
    !     |              |
    !     |              |
    !     ----------------
    !           (i,j)
end module Mesh

module ControlParameters
    use ConstantVariables
    use Mesh
    implicit none

    !--------------------------------------------------
    !Variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter                              :: CFL = 0.8 !CFL number
    ! real(KREAL), parameter                              :: MAX_TIME = 250.0 !Maximal simulation time
    integer(KINT), parameter                            :: MAX_ITER = 5E5 !Maximal iteration number
    real(KREAL), parameter                              :: EPS = 1.0E-5 !Convergence criteria
    real(KREAL)                                         :: simTime = 0.0 !Current simulation time
    integer(KINT)                                       :: iter = 1 !Number of iteration
    real(KREAL)                                         :: dt !Global time step
    real(KREAL)                                         :: res(4) !Residual
    
    !Output control
    character(len=6), parameter                         :: HSTFILENAME = "Cavity" !History file name
    character(len=6), parameter                         :: RSTFILENAME = "Cavity" !Result file name
    integer(KINT), parameter                            :: HSTFILE = 20 !History file ID
    integer(KINT), parameter                            :: RSTFILE = 21 !Result file ID

    !Gas propeties
    integer(KINT), parameter                            :: CK = 1 !Internal degree of freedom, here 1 denotes monatomic gas
    real(KREAL), parameter                              :: GAMMA = real(CK+4,KREAL)/real(CK+2,KREAL) !Ratio of specific heat
    ! real(KREAL), parameter                              :: OMEGA = 0.81 !Temperature dependence index in HS/VHS/VSS model
    real(KREAL), parameter                              :: PR = 2.0/3.0 !Prandtl number
    real(KREAL), parameter                              :: KN = 0.075 !Knudsen number in reference state
    real(KREAL), parameter                              :: ALPHA_REF = 1.0 !Coefficient in HS model
    real(KREAL), parameter                              :: OMEGA_REF = 0.5 !Coefficient in HS model
    real(KREAL), parameter                              :: MU_REF = 5.0*(ALPHA_REF+1.0)*(ALPHA_REF+2.0)*sqrt(PI)/(4.0*ALPHA_REF*(5.0-2.0*OMEGA_REF)*(7.0-2.0*OMEGA_REF))*KN !Viscosity coefficient in reference state

    !Geometry
    real(KREAL), parameter                              :: X_START = 0.0, X_END = 1.0, Y_START = 0.0, Y_END = 1.0 !Start point and end point in x, y direction 
    integer(KINT), parameter                            :: X_NUM = 45, Y_NUM = 45 !Points number in x, y direction
    integer(KINT), parameter                            :: IXMIN = 1 , IXMAX = X_NUM, IYMIN = 1 , IYMAX = Y_NUM !Cell index range
    integer(KINT), parameter                            :: N_GRID = (IXMAX-IXMIN+1)*(IYMAX-IYMIN+1) !Total number of cell
    integer(KINT), parameter                            :: GHOST_NUM = 1 !Ghost cell number

    !--------------------------------------------------
    !Initial flow field
    !--------------------------------------------------
    !Index method
    !-------------------------------
    !| (i-1) |  (i) |  (i) | (i+1) |
    !| cell  | face | cell | face  |
    !-------------------------------
    ! type(CellCenter)                                    :: ctr(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) !Cell center (with ghost cell)
    ! type(CellInterface)                                 :: vface(IXMIN-GHOST_NUM+1:IXMAX+GHOST_NUM),hface !Vertical and horizontal interfaces

    !Initial condition (density, u-velocity, v-velocity, lambda=1/temperature)
    real(KREAL), parameter, dimension(4)                :: INIT_GAS = [1.0, 0.0, 0.0, 1.0]

    !Boundary condition (density, u-velocity, v-velocity, lambda=1/temperature)
    real(KREAL), parameter, dimension(4)                :: BC_W = [1.0, 0.0, 0.0, 1.0] !West boundary
    real(KREAL), parameter, dimension(4)                :: BC_E = [1.0, 0.0, 0.0, 1.0] !East boundary
    real(KREAL), parameter, dimension(4)                :: BC_S = [1.0, 0.0, 0.0, 1.0] !South boundary
    real(KREAL), parameter, dimension(4)                :: BC_N = [1.0, 0.15, 0.0, 1.0] !North boundary

    !--------------------------------------------------
    !Discrete velocity space
    !--------------------------------------------------
    integer(KINT)                                       :: uNum = 64, vNum = 64 !Number of points in velocity space for u and v
    real(KREAL), parameter                              :: U_MIN = -15.0, U_MAX = +15.0, V_MIN = -15.0, V_MAX = +15.0 !Minimum and maximum micro velocity
    real(KREAL), allocatable, dimension(:,:)            :: uSpace,vSpace !Discrete velocity space for u and v
    real(KREAL), allocatable, dimension(:,:)            :: weight !Qudrature weight for discrete points in velocity space

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
        LocalFrame(3) = -w(2)*cosy+w(3)*cosx
        LocalFrame(4) = w(4)
    end function LocalFrame

    !--------------------------------------------------
    !>Convert macro variables from local frame to global
    !>@param[in] w            :macro variables in local frame
    !>@param[in] cosx,cosy    :directional cosine
    !>@return    GlobalFrame  :macro variables in global frame
    !--------------------------------------------------
    function GlobalFrame(w,cosx,cosy)
        real(KREAL) , intent(in)                        :: w(4)
        real(KREAL) , intent(in)                        :: cosx,cosy
        real(KREAL)                                     :: GlobalFrame(4)

        GlobalFrame(1) = w(1)
        GlobalFrame(2) = w(2)*cosx-w(3)*cosy
        GlobalFrame(3) = w(2)*cosy+w(3)*cosx
        GlobalFrame(4) = w(4)
    end function GlobalFrame

    !--------------------------------------------------
    !>Convert primary variables to conservative variables
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
    !>Convert conservative variables to primary variables
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
                    ((vn-prim(2))*qf(1)+(vt-prim(3))*qf(2))*(2*prim(4)*((vn-prim(2))**2+(vt-prim(3)**2))+CK-3)*B
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

        !allocate array
        allocate(sL(uNum,vNum))
        allocate(sR(uNum,vNum))

        sL = (midCell%h-leftCell%h)/(0.5*(midCell%length(idx)+leftCell%length(idx)))
        sR = (rightCell%h-midCell%h)/(0.5*(rightCell%length(idx)+midCell%length(idx)))
        midCell%sh(:,:,idx) = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

        sL = (midCell%b-leftCell%b)/(0.5*(midCell%length(idx)+leftCell%length(idx)))
        sR = (rightCell%b-midCell%b)/(0.5*(rightCell%length(idx)+midCell%length(idx)))
        midCell%sb(:,:,idx) = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)
    end subroutine VanLeerLimiter
end module Tools