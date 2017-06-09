# Tests #
Every test in __eiscor__ should adhere to the following template:
- a test may only directly test __one__ subroutine
- a name formatted as __test_prec_object_function__ (i.e. test_d_unifact_qr)
- a single call to the function __u_test_banner__ with the file name as an argument (i.e. __call u_test_banner('test_d_unifact_qr.f90')__)
- at every test for failure a call to the function __u_test_failed__ with a line number as an argument (i.e. __call u_test_failed(21)__)
- a single call at the end of the test to the function __u_test_passed__ with a total compute time as an argument (i.e. __call u_test_passed(0.0132)__)
- the last requirement implies that every test must also be timed
- if random tests are used, the seed must be explicitly set
- in DEBUG mode the test should check if the functions are responding correctly to incorrect inputs
- in VERBOSE mode the test may provide additional output

Below is a simple example of such a test that adheres to these requirements:
```fortran
!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_scalar_addition
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_scalar_addition

  implicit none
  
  ! compute variables
  real(8) :: a, b
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! set variables
  a = 2d0
  b = 3d0
    
  ! check +
  if ((a+b).NE.5d0) then
    call u_test_failed(__LINE__)
  end if
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_scalar_addition
```

## Executing tests ##
There are several convenient targets for building and running tests. The following targets require the installation of the __eiscor library__ and can be run from the __eiscor/__, __eiscor/test/__ and __eiscor/test/prec__ directories:
- __make tests__, compiles and executes all tests
- __make tests_d__, compiles and executes all double precision tests (does not work in __eiscor/test/complex_double/__)
- __make tests_z__, compiles and executes all complex double precision tests (does not work in __eiscor/test/double/__)
- __make tests_d_*__, compiles and executes all double precision tests containing __*__ (does not work in __eiscor/test/complex_double/__)
- __make tests_z_*__, compiles and executes all complex double precision tests __*__ (does not work in __eiscor/test/double/__)
