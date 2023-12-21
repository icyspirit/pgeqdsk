
      module fgeqdsk

      implicit none

      integer, parameter :: rkind = selected_real_kind(12,100)

      integer :: nw
      integer :: nh
      real(rkind) :: rdim
      real(rkind) :: zdim
      real(rkind) :: rcentr
      real(rkind) :: rleft
      real(rkind) :: zmid
      real(rkind) :: rmaxis
      real(rkind) :: zmaxis
      real(rkind) :: simag
      real(rkind) :: sibry
      real(rkind) :: bcentr
      real(rkind) :: current
      real(rkind), dimension(:), allocatable :: fpol
      real(rkind), dimension(:), allocatable :: pres
      real(rkind), dimension(:), allocatable :: ffprim
      real(rkind), dimension(:), allocatable :: pprime
      real(rkind), dimension(:,:), allocatable :: psirz
      real(rkind), dimension(:), allocatable :: qpsi
      integer :: nbbbs
      integer :: limitr
      real(rkind), dimension(:), allocatable :: rbbbs
      real(rkind), dimension(:), allocatable :: zbbbs
      real(rkind), dimension(:), allocatable :: rlim
      real(rkind), dimension(:), allocatable :: zlim

      contains

      subroutine read_geqdsk(filename)

      implicit none

      character(*) :: filename

      integer :: idum
      real(rkind) :: xdum
      character*10 zcase(6)
      integer :: neqdsk,i,j

      open(neqdsk, file=filename)

      read (neqdsk,2000) (zcase(i),i=1,6),idum,nw,nh
      read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
      read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
      read (neqdsk,2020) current,simag,xdum,rmaxis,xdum
      read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      if (allocated(fpol)) deallocate(fpol)
      if (allocated(pres)) deallocate(pres)
      if (allocated(ffprim)) deallocate(ffprim)
      if (allocated(pprime)) deallocate(pprime)
      if (allocated(psirz)) deallocate(psirz)
      if (allocated(qpsi)) deallocate(qpsi)
      allocate(fpol(nw))
      allocate(pres(nw))
      allocate(ffprim(nw))
      allocate(pprime(nw))
      allocate(psirz(nw,nh))
      allocate(qpsi(nw))
      read (neqdsk,2020) (fpol(i),i=1,nw)
      read (neqdsk,2020) (pres(i),i=1,nw)
      read (neqdsk,2020) (ffprim(i),i=1,nw)
      read (neqdsk,2020) (pprime(i),i=1,nw)
      read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
      read (neqdsk,2020) (qpsi(i),i=1,nw)
      read (neqdsk,2022) nbbbs,limitr
      if (allocated(rbbbs)) deallocate(rbbbs)
      if (allocated(zbbbs)) deallocate(zbbbs)
      if (allocated(rlim)) deallocate(rlim)
      if (allocated(zlim)) deallocate(zlim)
      allocate(rbbbs(nbbbs))
      allocate(zbbbs(nbbbs))
      allocate(rlim(limitr))
      allocate(zlim(limitr))
      read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)

      close(neqdsk)

 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)

      end subroutine read_geqdsk

      subroutine write_geqdsk(filename)

      implicit none

      character(*) :: filename

      integer :: idum = 0
      real(rkind) :: xdum = 0.0
      integer :: neqdsk,i,j

      if (.not.allocated(fpol)) call abort
      if (.not.allocated(pres)) call abort
      if (.not.allocated(ffprim)) call abort
      if (.not.allocated(pprime)) call abort
      if (.not.allocated(psirz)) call abort
      if (.not.allocated(qpsi)) call abort
      if (.not.allocated(rbbbs)) call abort
      if (.not.allocated(zbbbs)) call abort
      if (.not.allocated(rlim)) call abort
      if (.not.allocated(zlim)) call abort

      open(neqdsk, file=filename, status='unknown')

      write (neqdsk,2000) ('GEQDSK',i=1,6),idum,nw,nh
      write (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
      write (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
      write (neqdsk,2020) current,simag,xdum,rmaxis,xdum
      write (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      write (neqdsk,2020) (fpol(i),i=1,nw)
      write (neqdsk,2020) (pres(i),i=1,nw)
      write (neqdsk,2020) (ffprim(i),i=1,nw)
      write (neqdsk,2020) (pprime(i),i=1,nw)
      write (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
      write (neqdsk,2020) (qpsi(i),i=1,nw)
      write (neqdsk,2022) nbbbs,limitr
      write (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      write (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)

      close(neqdsk)

 2000 format(6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)

      end subroutine write_geqdsk

      end module
