!> https://fortran-lang.discourse.group/t/best-way-to-declare-a-double-precision-in-fortran/69/2
module mod_kinds
use iso_fortran_env, only:  real32, real64
implicit none
integer, parameter :: sp = real32
integer, parameter :: dp = real64
end module mod_kinds
