! Shared constants

module ConstantsModule

  implicit none

  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)
  real(sp), parameter :: pi = 3.141592653589793238462643383279502884197_sp
  real(dp), parameter :: pi_d = 3.141592653589793238462643383279502884197_dp

end module ConstantsModule
