program tutorial
    implicit none
    real :: rand(1)
    call random_number(rand)
    print "(a20,i1)","ABS(-1) = ", ABS(-1)
    print "(a20,i1)","INT(2.71) = ", INT(2.71)
    print "(a20,i1)","NINT(2.71) = ", NINT(2.71)
    print "(a20,i1)","FLOOR(2.71) = ", FLOOR(2.71)
    print "(a20,f6.3)","COS(0.0) = ", COS(0.0)
    print "(a20,f6.3)","COS(1.0) = ", COS(1.0)
    print "(a20,f6.3)","COS(2.0) = ", COS(2.0)
    print "(a20,f6.3)","COS(3.14159) = ", COS(3.14159)
    print "(a20,f6.3)","COS(2*pi) = ", COS(2*3.14159)


end program tutorial
