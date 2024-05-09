subroutine simple(a,b,c)

    real, intent(in) :: a, b
    real, intent(out) :: c

    c = a + b

end subroutine simple

function foo(a) result(b)
    implicit none
    
    real(kind=8), intent(in)    :: a(:,:)
    complex(kind=8)             :: b(size(a,1),size(a,2))
    
    b = exp((0,1)*a)
    
end function foo