module header_mod

contains

subroutine pine_tree()

    implicit none
    character(len=70), dimension(13) :: tree = [ &
        "⠀⠀⠀⢀⣀⡴⠲⡄⣀⠀⣠⠖⢆⣀⡤⣄⣠⣄⠀⠀⠀⠀⠈", &
        "⠀⠀⡤⠷⠄⠀⠀⠀⠀⠙⠋⠀⠀⠀⠀⠀⠈⠁⠀⠐⡄⠀⠀", &
        "⠀⡜⠀⠀⠤⠊⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠄⢤⠀⠘⡄", &
        "⢠⣇⠀⠀⠀⠤⠶⠶⡒⣋⣯⠖⡇⣀⡀⠀⠀⠀⠀⠀⡀⡠⠇", &
        "⠀⠀⠀⠉⠉⠲⠶⢖⡟⠋⢀⡔⠃⡏⠦⠄⠄⠤⠴⠃⠀⠀⠀", &
        "⠀⠀⠀⠀⠀⠀⠀⠀⡧⠖⠋⠀⡸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀", &
        "⠀⠀⠀⠀⠀⠀⠀⡞⠀⠀⡠⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀", &
        "⠀⠀⠀⠀⠀⠀⢸⠀⠀⠘⠲⠤⣠⠶⡄⠀⠀⠀⠀⠀⠀⠀⠀", &
        "⠀⠀⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⠀⡏⠉⠀⠀⠀⠀⠀⠀⠀⠀", &
        "⠀⠀⠀⠀⠀⠑⠢⣄⡀⠀⠀⠀⠀⡇⠀⠀⠈⠉⠉⠁⠀⠀⠀", &
        "⠀⠀⠀⠀⠀⠀⠀⠀⠈⡇⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀", &
        "⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⣸⠀⠀⠀⠐⠋⠀⠀⠀⠀⠀", &
        "⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀" ]

    integer :: i
    do i = 1, size(tree)
        write(*,*) trim(tree(i))
    end do
    write(*,*) ""
    write(*,*) ""
end subroutine 

subroutine final_phrase()
    
    implicit none
    character(len=:), allocatable :: phrases(:)
    real :: r
    integer :: idx
    character(len=256) :: output

    phrases = [character(len=256) :: &
    "I'll break it down just to feel okay, I'll tear it down just to feel the pain."&
    & //achar(10)// "Sorrower, ~ Orbit Culture", &
    
    "I'm scared to get close, and I hate being alone." & 
    &//achar(10)// "Can You Feel My Heart?, Bring Me The Horizon", &
    "Everybody showed up for the execution but nobody would show their face." & 
    &//achar(10)// "One Assasination Under God, Marilyn Manson", &
    
    "The only dream that matters is the one you wake up from."& 
    &//achar(10)// "Faithless by Default, Dark Tranquillity", &
    
    "Fear is the weakness in all of us, it's sad to see you go."& 
    &//achar(10)// "Fear Is The Weakness, In Flames",&
    
    "If this is my crown, I am its slave." & 
    &//achar(10)// "Death Can Take Me, Lorna Shore",&

    "Wenn Fliegen hinter Fliegen fliegen, fliegen Fliegen hinter Fliegen her."&
    &//achar(10)// "~A German tongue twister",&

    "Thanks to quantum mechanics, we seem to do have a free will, but I prefer the deterministic world"&
    &//achar(10)// "~based on a personal discussion (but Einstein and Bohr had that discussion too)",&

    "Pine trees are the best trees, that is a fact."&
    &//achar(10)// "~(definitely not my subjective opinion)",&


    "We're dancing like flames flickering in the night"& 
    &//achar(10)// "Pain Remains I, Lorna Shore" ]


    ! Initialize the random number generator
    call random_seed()

    ! Generate a random number 0 ≤ r < 1
    call random_number(r)

    ! Convert to an index from 1 to size(phrases)
    idx = int(r * size(phrases)) + 1
    output = char(27)//"[3m"//trim(phrases(idx))//char(27)//"[23m"
    write(*,*) ""
    write(*,*) ""
    write(*,*) output
end subroutine 


end module
