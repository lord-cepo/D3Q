module delta_nodes_bind
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine c_delta_nodes(values, gradients, quality, qorder, number_of_nodes, weights, coords) &
      bind(C, name="delta_nodes")
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: values
      type(c_ptr), value :: gradients
      integer(c_int), value :: quality, qorder
      integer(c_int), intent(out) :: number_of_nodes
      type(c_ptr), intent(out) :: weights
      type(c_ptr), intent(out) :: coords
    end subroutine c_delta_nodes
  end interface

  Interface

  Subroutine c_free(ptr) BIND(C,name="free_c_ptr")
    IMPORT :: C_PTR
    Implicit NONE
    Type(C_PTR) :: ptr
  End Subroutine c_free

End Interface
CONTAINS
  subroutine delta_nodes_f(values, gradients, quality, qorder, number_of_nodes, weights, coords)
    use kinds, only: dp
    use iso_c_binding
    real(dp), INTENT(IN) :: values(8), gradients(3,8)
    integer, INTENT(IN) :: quality, qorder
    integer, INTENT(OUT) :: number_of_nodes
    real(dp), ALLOCATABLE, INTENT(OUT) :: weights(:), coords(:,:)

    integer :: i,j
    integer(c_int) :: quality_, qorder_, number_of_nodes_
    real(c_double), dimension(:), pointer :: values_
    real(c_double), dimension(:), pointer :: gradients_
    real(c_double), dimension(:), pointer :: weights_
    real(c_double), dimension(:), pointer :: coords_
    type(c_ptr) :: c_values, c_gradients, c_weights, c_coords

    allocate(values_(8))
    allocate(gradients_(24))
    values_ = REAL(values, kind=c_double)
    FORALL(j=1:8, i=1:3) gradients_(3*j + i + 1) = REAL(gradients(i,j), kind=c_double)
    quality_ = INT(quality, kind=c_int)
    qorder_  = INT(qorder,  kind=c_int)

    c_values = c_loc(values_)
    c_gradients = c_loc(gradients_)

    call c_delta_nodes(c_values, c_gradients, quality_, qorder_, number_of_nodes_, c_weights, c_coords)
    number_of_nodes = INT(number_of_nodes_)
    call c_f_pointer(c_weights, weights_, [number_of_nodes])
    call c_f_pointer(c_coords, coords_, [3*number_of_nodes])

    ! Clean up
    allocate(weights(number_of_nodes))
    weights = REAL(weights_, kind=dp)
    allocate(coords(3,number_of_nodes))
    coords = REAL(RESHAPE(coords_, SHAPE(coords)), kind=dp)
    
    deallocate(values_)
    deallocate(gradients_)
    
    ! deallocate(weights_)
    ! deallocate(coords_)
    ! call c_free(c_values)
    ! call c_free(c_gradients)
    ! call c_free(c_weights)
    ! call c_free(c_coords)

  end subroutine delta_nodes_f
end module