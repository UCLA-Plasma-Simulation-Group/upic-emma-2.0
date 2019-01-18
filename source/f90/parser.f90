! ! Written by Han Wen, UCLA
!
! A recursive descent parser for the following parsing expression grammar (PEG)
!
! Opexpr  = L4EXPR EOL
! L0EXPR  = L1EXPR (L1OPR L1EXPR)*
! L1EXPR  = L2EXPR (L2OPR L2EXPR)*
! L2EXPR  = L3OPR? L3EXPR
! L3EXPR  = L4EXPR (L4OPR L4EXPR)?
! L4EXPR  = L4OPR? L5EXPR (L5OPR L5EXPR)*
! L5EXPR  = L6EXPR (L6OPR L6EXPR)*
! L6EXPR  = term (L7OPR L6EXPR)*
! term    = func / brt / var / num
! brt     = '(' L0EXPR ')'
! func    = ~"[a-z A-Z _]+(([0-9]*)([a-z A-Z _]*))*" '(' L0EXPR (',' L0EXPR)* ')'
! var     = ~"[a-z A-Z _]+[0-9]+"
! num     = ~"[\-+]?([0-9]*)?(\.)?([0-9]*)?([dDeE][\-+]?[0-9]+)?"
! L1OPR   = '||'
! L2OPR   = '&&'
! L3OPR   = '~' / '!'
! L4OPR   = '<=' / '<' / '>=' / '>' / '==' / '!=' / '/='
! L5OPR   = '+' / '-'
! L6OPR   = '*' / '/'
! L7OPR   = '**' / '^'
! 

! basic ideas: string => AST => register based bytecode (easier CSE?) => VM execution
! compile time is not critical but we want fast execution. So we need some
! optimizations when compiling (for example, common subexpression elimination, branch
! elimination, constant folding. Even maybe expression simplication)

! this module can be tested alone by uncomenting the main routine at the end of this file

module m_rvm
    implicit none
    private
    integer, parameter, public  :: p_k_parse = selected_real_kind(15)
    integer, parameter, public  :: MAX_LEN_PSR = 2048  ! maximum length for the input string
    integer, parameter          :: MAX_LEN_OPR = 16  ! maximum length of the name of an instructions (function)
    integer, parameter          :: MAX_LEN_VAR = 8  ! maximum length of the name of variable
    integer, parameter          :: MAX_NUM_OPRD = 4  ! maximum number of operands in one instruction
    integer, parameter          :: NUM_FUNCTIONS= 31  ! number of buildin functions
    integer, parameter          :: NUM_OPRS = 16  ! number of buildin operators
    integer, parameter          :: NUM_RULES = 8  ! number of nonterminal rules
    integer, parameter          :: NUM_CONSTANT = 4 ! number of contants
    logical, parameter          :: DEBUG = .true.  ! enable debuging
    logical, parameter          :: DEBUG_VERBOSE = .false.  ! more debug info
    logical, parameter          :: DEBUG_RVM = .true.  ! VM status for human reader

    !---------------------------------------------------------------------------------------------------
    !===================================================================================================
    !---------------------------------------------------------------------------------------------------
    type :: t_buildin
        character (len=MAX_LEN_OPR) :: name = ''
        integer                     :: n_operands = 0
        integer                     :: id = 0
    end type t_buildin
 
    type :: t_constant
        character (len=MAX_LEN_OPR) :: name = ''
        real(p_k_parse)             :: val = 0.
    end type t_constant

    integer, parameter  :: ILEN = MAX_NUM_OPRD + 1
    type :: t_parse_result
        ! instructions format: /op_id, destination, source1, source2 .../ where
        ! op_id is the id of the instruction
        ! destination is the register address where the result be written to
        ! source* are the register address from which to read operands/parameters
        ! note that destination/source* /= 0, and 0 suggests the end of the instruction
        integer, dimension(ILEN) :: code = 0
        
        integer                  :: next = 1  ! next position to parse
    end type t_parse_result

    type :: t_parser
        character (len=MAX_LEN_PSR) :: expr = ""
        integer                     :: slen = 0  ! length of the string not counting whitespaces
        integer                     :: pp   = 1  ! current parsing position
        integer                     :: rp   = 1  ! current register position
        integer                     :: ip   = 0  ! current instruction queue position
        integer                     :: nvar = 0  ! number of variables
        character (len=MAX_LEN_VAR), dimension(:), allocatable :: var_name
        integer, dimension(:, :), allocatable :: inst
        integer, dimension(:, :), allocatable :: inst_parse  ! store register value at parse (tmp)
        real(p_k_parse), dimension(:), allocatable :: reg
        real(p_k_parse), dimension(:), allocatable :: reg_parse  ! store register value at parse (tmp)
        ! reversely linked to the instruction in memo writting to this register
        !integer, dimension(:,:), allocatable :: reg_parse_link
    end type t_parser

    ! compile time constants are .le. 0
    integer, parameter :: flag_noparse  = 0,  &
                          flag_load     = -1, &
                          flag_success  = -2, &
                          flag_failed   = -3, &
                          p_error       = -100, &
                          err_num       = p_error-1, &
                          err_func_nf   = p_error-2, &
                          err_func_ipc  = p_error-3, &
                          err_brt_mtch  = p_error-4, &
                          err_brt_cntt  = p_error-5, &
                          err_var       = p_error-6, &
                          err_opr       = p_error-7, &
                          err_eof       = p_error-8, &
                          err_noterm    = p_error-9, &
                          err_runtime   = p_error-9999

    ! runtime instructions
    ! NOTE1: if you want to add a function, you need to modify 4 places:
    ! 1) increase the NUM_FUNCTIONS parameter above, 2) add the signature below, 
    ! 3) add signature to the functions list below, 4) add definition to the eval_parse function
    ! 
    ! NOTE2: if you want to add a constant, you need to modify 2 places:
    ! 1) increase the NUM_CONSTANT parameter above
    ! 2) add constants to the list below, note that it has lower priority than user
    ! defined variables
    type(t_buildin), parameter :: op_noop      = t_buildin(''             , 0, 1),  &
                                  op_or        = t_buildin('||'           , 2, 2),  &
                                  op_and       = t_buildin('&'//'&'       , 2, 3),  & ! &&
                                  op_not       = t_buildin('! ~'          , 2, 4),  &
                                  op_le        = t_buildin('<='           , 2, 5),  &
                                  op_lt        = t_buildin('<'            , 2, 6),  &
                                  op_ge        = t_buildin('>='           , 2, 7),  &
                                  op_gt        = t_buildin('>'            , 2, 8),  &
                                  op_eq        = t_buildin('=='           , 2, 9),  &
                                  op_ne        = t_buildin('!= /='        , 2, 10),  &
                                  op_plus      = t_buildin('+'            , 2, 11), &
                                  op_minus     = t_buildin('-'            , 2, 12), &
                                  op_mul       = t_buildin('*'            , 2, 13), &
                                  op_div       = t_buildin('/'            , 2, 14), &
                                  op_mod       = t_buildin('%'            , 2, 15), &
                                  op_pow       = t_buildin('^ **'         , 2, 16), &
                                  func_abs     = t_buildin('abs'          , 1, 17) , &
                                  func_sinh    = t_buildin('sinh'         , 1, 18) , &
                                  func_sin     = t_buildin('sin'          , 1, 19) , &
                                  func_cosh    = t_buildin('cosh'         , 1, 20) , &
                                  func_cos     = t_buildin('cos'          , 1, 21) , &
                                  func_tanh    = t_buildin('tanh'         , 1, 22) , &
                                  func_tan     = t_buildin('tan'          , 1, 23) , &
                                  func_exp     = t_buildin('exp'          , 1, 24) , &
                                  func_log10   = t_buildin('log10'        , 1, 25) , &
                                  func_log     = t_buildin('log'          , 1, 26) , &
                                  func_asin    = t_buildin('asin arcsin'  , 1, 27) , &
                                  func_acos    = t_buildin('acos arccos'  , 1, 28) , &
                                  func_atan2   = t_buildin('atan2 arctan2', 2, 29) , &
                                  func_atan    = t_buildin('atan arctan'  , 1, 30) , &
                                  func_sqrt    = t_buildin('sqrt'         , 1, 31) , &
                                  func_not     = t_buildin('not'          , 1, op_not%id) , &
                                  func_neg     = t_buildin('neg'          , 1, 33) , &
                                  func_if      = t_buildin('if'           , 3, 34) , &
                                  func_pow     = t_buildin('pow'          , 2, 35) , &
                                  func_int     = t_buildin('int'          , 1, 36) , &
                                  func_nint    = t_buildin('nint'         , 1, 37) , &
                                  func_ceiling = t_buildin('ceiling'      , 1, 38) , &
                                  func_floor   = t_buildin('floor'        , 1, 39) , &
                                  func_modulo  = t_buildin('modulo'       , 2, 40) , &
                                  func_mod     = t_buildin('mod'          , 2, op_mod%id) , &
                                  func_rect    = t_buildin('rect'         , 1, 41) , &
                                  func_step    = t_buildin('step'         , 1, 42) , &
                                  func_min3    = t_buildin('min3'         , 3, 43) , &
                                  func_min     = t_buildin('min'          , 2, 44) , &
                                  func_max3    = t_buildin('max3'         , 3, 45) , &
                                  func_max     = t_buildin('max'          , 2, 46)  

    type(t_buildin), dimension(NUM_FUNCTIONS), parameter :: functions = (/ func_abs, func_sinh, func_sin, func_cosh, func_cos, &
            func_tanh, func_tan, func_exp, func_log10, func_log, func_asin, func_acos, func_atan2, func_atan, func_sqrt, &
            func_not, func_neg, func_if, func_pow, func_int, func_nint, func_ceiling, func_floor, func_modulo, func_mod, &
            func_rect, func_step, func_min3, func_min, func_max3, func_max/) 
                                                        
    type(t_constant), dimension(NUM_CONSTANT ), parameter :: constants = (/ &
                                  t_constant('pi           ', 3.141592653589793238462643383_p_k_parse ), &
                                  t_constant('e            ', 2.718281828459045235360287471_p_k_parse ), &
                                  t_constant('true         ', 1.0_p_k_parse                           ), &
                                  t_constant('false        ', 0.0_p_k_parse                           ) /)

   type(t_buildin), dimension(NUM_OPRS), parameter :: operators = (/ op_noop, op_or, op_and, op_not, op_le, op_lt, &
       op_ge, op_gt, op_eq, op_ne, op_plus, op_minus, op_mul, op_div, op_mod, op_pow /)






    interface setup
        module procedure setup_parser
    end interface
    
    interface delete
        module procedure delete_parser
    end interface

    interface eval
        module procedure eval_parser
    end interface
    
    interface functext
        module procedure functext_parser
    end interface

    public :: t_parser, setup, eval, functext, delete

    contains

!----------------------------------------------------------------------------------------------------------
! return the string parser sees
!----------------------------------------------------------------------------------------------------------
    function functext_parser (this)
        implicit none
        type(t_parser),  intent(in) :: this
        character (len=len_trim(this%expr)) :: functext_parser

        functext_parser = trim(this%expr)
    end function

!----------------------------------------------------------------------------------------------------------
! free resources after parsing
!----------------------------------------------------------------------------------------------------------
    subroutine cleanup_parser (this)
        implicit none
        type(t_parser),  intent(inout) :: this

        !deallocate ( this%reg_parse_link )
        deallocate ( this%reg_parse )
        deallocate ( this%inst_parse )
    end subroutine

!----------------------------------------------------------------------------------------------------------
! free all memory 
!----------------------------------------------------------------------------------------------------------
    subroutine delete_parser (this)
        implicit none
        type(t_parser),  intent(inout) :: this

        deallocate ( this%reg )
        deallocate ( this%inst )
        deallocate ( this%var_name )
    end subroutine

!----------------------------------------------------------------------------------------------------------
! test if this%expr contains a buildin instruction at current position, 
! return the length of the matched operator string; return 0 if not
!----------------------------------------------------------------------------------------------------------
    function is_buildin( this, op )
        implicit none
        type(t_parser),  intent(in) :: this
        type(t_buildin), intent(in) :: op
        integer :: is_buildin

        !local
        integer :: i, j, e
        is_buildin = 0
        e = 1
        do while (e <= len_trim(op%name))
            j = e
            ! find next white space
            do i = e, len_trim(op%name)+1
                if (op%name(i:i) == ' ') then
                    e = i - 1
                    exit
                endif
            enddo
            !if (DEBUG_VERBOSE) print*, 'try to match ', this%expr(this%pp:this%pp+e-j), ' with ', op%name(j:e)
            if (op%name(j:e) == this%expr(this%pp:this%pp+e-j)) then
                is_buildin = e - j + 1
                if (DEBUG) print*, 'matched opr "', op%name(j:e), '", next pos = ', this%pp+e-j+1
                exit
            endif
            ! did not find one, proceed to next pattern
            e = e + 2
        enddo
    end function
    
!----------------------------------------------------------------------------------------------------------
! record the parsing result
!----------------------------------------------------------------------------------------------------------
    subroutine write_stat( this, stat, next, co )
        implicit none
        type(t_parser),        intent (inout) :: this
        type(t_parse_result),  intent (inout) :: stat
        integer,               intent (in)    :: next ! next: next position to parse
        integer, dimension(:), intent (in)    :: co  ! instruction to record
        stat%next = next
        stat%code(1:size(co)) = co
        if (co(1) > flag_failed) then
            !this%reg_parse_link(:, this%rp) = (/ level, this%pp /)
            if (co(1) /= op_noop%id) then
                ! if loading the variables or constants, don't increase register count
                if (.not. (co(1) == flag_load .and. co(2) <= this%nvar + NUM_CONSTANT)) &
                    & this%rp = this%rp + 1
                this%ip = this%ip + 1
                this%inst_parse(:, this%ip) = stat%code
            endif
        endif
    end subroutine

!----------------------------------------------------------------------------------------------------------
! match the pattern: PATTERN1 (OPR PATTERN2)*
!----------------------------------------------------------------------------------------------------------
    recursive subroutine meta_rule1( this, pattern1, pattern2, opr, stat, times )
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        type(t_buildin), dimension(:),     intent (in)    :: opr
        integer, optional,    intent (in)    :: times
        interface

            recursive subroutine pattern1(this, stat)
                import t_parser
                import t_parse_result
                type(t_parser),       intent (inout) :: this
                type(t_parse_result), intent (inout) :: stat
            end subroutine

            recursive subroutine pattern2(this, stat)
                import t_parser
                import t_parse_result
                type(t_parser),       intent (inout) :: this
                type(t_parse_result), intent (inout) :: stat
            end subroutine
        end interface
        
        ! local
        integer :: pp, i, t, maxt, d
        type(t_parse_result) :: lstat, rstat
        if (present(times)) then
            maxt = times
        else
            maxt = len_trim(this%expr)
        endif
        pp = this%pp  ! back up the position because this%pp may be changed by term ()

        call pattern1(this, lstat)
        if (lstat%code(1) > flag_failed) then  ! left operand is legit
            do t = 1, maxt
                ! find the operator
                d = 0
                do i = 1, size(opr)
                    d = is_buildin( this, opr(i))
                    if ( d > 0 ) then
                        this%pp = this%pp + d  ! consume operator
                        exit
                    endif
                enddo
                if ( d > 0) then
                    call pattern2(this, rstat)
                    if (rstat%code(1) == flag_failed) then  ! second pattern not found
                        !call handle_error(this, this%memo(level, this%pp)%code, (/ err_opr, this%pp /), trim(opr(i)%name) )
                        call handle_error(this, stat%code, (/ err_opr, this%pp /), trim(opr(i)%name) )
                    else                                ! we match the pattern successfully
                        call write_stat( this, stat, rstat%next, (/ opr(i)%id, this%rp, lstat%code(2), rstat%code(2) /))
                        if (DEBUG) print*, 'matched "', this%expr(pp:this%pp-1), '", instruction =', stat%code(1:4), stat%next
                        lstat = stat
                    endif
                else
                    exit
                endif
            enddo
        endif
        stat = lstat
    end subroutine meta_rule1

!----------------------------------------------------------------------------------------------------------
! match pattern OPR? PATTERN
!----------------------------------------------------------------------------------------------------------
    recursive subroutine meta_rule2( this, pattern, opr, stat, opr_reg )
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        type(t_buildin), dimension(:),     intent (in)    :: opr
        type(t_buildin), dimension(:),optional, intent (in)    :: opr_reg
        interface
            recursive subroutine pattern(this, stat)
                import t_parser
                import t_parse_result
                type(t_parser),       intent (inout) :: this
                type(t_parse_result), intent (inout) :: stat
            end subroutine
        end interface
        
        !local
        integer :: d, i, pp
        d = 0
        pp = this%pp
        do i = 1, size(opr)
            d = is_buildin( this, opr(i))
            if ( d > 0 ) then
                this%pp = this%pp + d  ! consume operator
                exit
            endif
        enddo
        call pattern(this, stat)
        if ( d > 0 ) then
            if (stat%code(1) == flag_failed) then
                !call handle_error( this, this%memo(level, this%pp)%code, (/ err_opr, pp /), trim(opr(i)%name) )
                call handle_error( this, stat%code, (/ err_opr, pp /), trim(opr(i)%name) )
            else
                if (present(opr_reg)) then
                    call write_stat( this, stat, stat%next, (/ opr_reg(i)%id, this%rp, stat%code(2) /))
                else
                    call write_stat( this, stat, stat%next, (/ opr(i)%id, this%rp, stat%code(2) /))
                endif
                if (DEBUG) print*, 'matched "', this%expr(pp:this%pp-1), '", instruction =', stat%code(1:3), stat%next
            endif
        endif
    end subroutine

!----------------------------------------------------------------------------------------------------------
! L6EXPR = term ([^ **] L6EXPR)*
!----------------------------------------------------------------------------------------------------------
    recursive subroutine l6expr(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        call meta_rule1(this, term, l6expr, (/op_pow/), stat)
        if (DEBUG_VERBOSE) print*, 'l6expr stat=', stat
    end subroutine

!----------------------------------------------------------------------------------------------------------
! L5EXPR = L6EXPR ([*/%] L6EXPR)*
!----------------------------------------------------------------------------------------------------------
    recursive subroutine l5expr(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        call meta_rule1(this, l6expr, l6expr, (/op_mul, op_div, op_mod/), stat)
        if (DEBUG_VERBOSE) print*, 'l5expr stat=', stat
    end subroutine

!----------------------------------------------------------------------------------------------------------
! match [-+]? L5EXPR
!----------------------------------------------------------------------------------------------------------
    recursive subroutine np_l5expr(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        call meta_rule2(this, l5expr, (/op_plus, op_minus/), stat, (/op_noop, func_neg/))
        if (DEBUG_VERBOSE) print*, 'np_l5expr stat=', stat
    end subroutine

!----------------------------------------------------------------------------------------------------------
! L4EXPR  = [-+]? L5EXPR ([-+] L5EXPR)*
!----------------------------------------------------------------------------------------------------------
    recursive subroutine l4expr(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        call meta_rule1(this, np_l5expr, l5expr, (/op_plus, op_minus/), stat)
        if (DEBUG_VERBOSE) print*, 'l4expr stat=', stat
    end subroutine

!----------------------------------------------------------------------------------------------------------
! L3EXPR  = L4EXPR ([<= < >= > == != /=] L4EXPR)?
! here we forbid something like 0<x<=1 and x>0 == y<0, use 0<x && x<=1 and (x>0) == (y<0) instead
!----------------------------------------------------------------------------------------------------------
    recursive subroutine l3expr(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        call meta_rule1(this, l4expr, l4expr, (/op_le, op_lt, op_ge, op_gt, op_eq, op_ne/), stat, 1)
        if (DEBUG_VERBOSE) print*, 'l3expr stat=', stat
    end subroutine

!----------------------------------------------------------------------------------------------------------
! L2EXPR  = [! ~]? L3EXPR
!----------------------------------------------------------------------------------------------------------
    recursive subroutine l2expr(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        call meta_rule2(this, l3expr, (/op_not/), stat)
        if (DEBUG_VERBOSE) print*, 'l2expr stat=', stat
    end subroutine

!----------------------------------------------------------------------------------------------------------
! L1EXPR  = L2EXPR ('&&' L2EXPR)*
!----------------------------------------------------------------------------------------------------------
    recursive subroutine l1expr(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        call meta_rule1(this, l2expr, l2expr, (/op_and/), stat)
        if (DEBUG_VERBOSE) print*, 'l1expr stat=', stat
    end subroutine

!----------------------------------------------------------------------------------------------------------
! L0EXPR  = L1EXPR (L1OPR L1EXPR)*
!----------------------------------------------------------------------------------------------------------
    recursive subroutine l0expr(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        call meta_rule1(this, l1expr, l1expr, (/op_or/), stat)
        if (DEBUG_VERBOSE) print*, 'l0expr stat=', stat
    end subroutine

!----------------------------------------------------------------------------------------------------------
! match the user defined variable name
!----------------------------------------------------------------------------------------------------------
    subroutine variable(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        
        ! local
        character :: c
        integer :: pp, i
        logical :: found
        found = .false.
        do i = 1, size(this%var_name)
            pp = this%pp + len_trim(this%var_name(i)) - 1
            if (trim(this%var_name(i)) == this%expr(this%pp:pp)) then
                found = .true.
                exit
            endif
        enddo
        if (found) then
            call write_stat( this, stat, pp+1, (/ flag_load, i /))
            this%pp = pp + 1
            if (DEBUG) print*, 'matched variable "', trim(this%var_name(i)), '", next pos = ', this%pp
        else
            call write_stat( this, stat, this%pp, (/ flag_failed /))
        endif
    end subroutine

!----------------------------------------------------------------------------------------------------------
! match numbers [\-+]?([0-9]*)?(\.)?([0-9]*)?([dDeE][\-+]?[0-9]+)?
!----------------------------------------------------------------------------------------------------------
    subroutine number(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat
        
        ! local
        character :: c
        integer :: pp
        logical :: deci, expo
        real(p_k_parse) :: val
        deci = .False.
        expo = .False.
        pp = this%pp  ! get original position
        c = this%expr(pp:pp)
        if (scan(c, '+-0123456789.') > 0) then
            if (c == '+' .or. c == '-') pp = pp + 1
            do while (pp <= this%slen)
                if (c == '.') then
                    if (deci .or. expo) then
                        call handle_error(this, stat%code, (/err_num, pp/))
                    else
                        deci = .true.
                        if ( '0' > this%expr(pp-1:pp-1) .or. &
                                this%expr(pp-1:pp-1) > '9') then 
                            if ( '0' > this%expr(pp+1:pp+1) .or. &
                                    this%expr(pp+1:pp+1) > '9') then 
                                ! no digit around '.', raise error
                                call handle_error(this, stat%code, (/err_num, pp/))
                            endif
                            pp = pp + 1
                        endif
                    endif
                else if (c == 'e' .or. c == 'd') then
                    if (expo .eqv. .false.) then
                        expo = .true.
                        ! look at next one or position
                        if (scan(this%expr(pp+1:pp+1), '+-') > 0) then
                            if (scan(this%expr(pp+2:pp+2), '0123456789') > 0) then
                                pp = pp + 2
                            else
                                call handle_error(this, stat%code, (/err_num, pp+2/))
                            endif
                        else if (scan(this%expr(pp+1:pp+1), '0123456789') > 0) then
                            pp = pp + 1
                        else
                            call handle_error(this, stat%code, (/err_num, pp+1/))
                        endif
                    else
                        call handle_error(this, stat%code, (/err_num, pp/))
                    endif
                else
                    if ( '0' > c .or. c > '9') exit
                endif
                pp = pp + 1
                c = this%expr(pp:pp)
            enddo
            if (DEBUG) print*, 'matched number "', this%expr(this%pp:pp-1), '", save to register', this%rp
            !---------------- mixing parsing and compiling codes ------------------
            read(this%expr(this%pp:pp-1), *) val
            this%reg_parse(this%rp) = real(val, p_k_parse)
            !---------------- mixing parsing and compiling codes ------------------
            call write_stat( this, stat, pp, (/ flag_load, this%rp, this%pp /))
            this%pp = pp  ! update current pos
        else  ! not a number
            call write_stat( this, stat, this%pp, (/ flag_failed /))
        endif
    end subroutine

!----------------------------------------------------------------------------------------------------------
! term    =  func / num / bracket / var / constant
!----------------------------------------------------------------------------------------------------------
    recursive subroutine term( this, stat )
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat

        ! early termination if we are at the end of the string
        if (len_trim(this%expr) < this%pp) then
            stat%code(1) = flag_failed
            return
        endif

        ! check out term one by one
        call func(this, stat)

        if (stat%code(1) == flag_failed) call number(this, stat)

        if (stat%code(1) == flag_failed) call bracket(this, stat)

        if (stat%code(1) == flag_failed) call variable(this, stat)

        if (stat%code(1) == flag_failed) call constant(this, stat)

        if (stat%code(1) == flag_failed) call handle_error(this, stat%code, (/err_noterm, this%pp/))
    end subroutine


!----------------------------------------------------------------------------------------------------------
! func    = ~"[a-z A-Z _]+(([0-9]*)([a-z A-Z _]*))*" '(' L0EXPR (',' L0EXPR)* ')'
!----------------------------------------------------------------------------------------------------------
    recursive subroutine func(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat

        ! local
        integer :: j, np, pp, d
        integer, dimension(MAX_NUM_OPRD) :: mo
        character :: c
        pp = this%pp
        ! ! this is a more general match
        !c = this%expr(this%pp:this%pp)
        !if( ( c >= 'a' .and. c <= 'z') .or. c == '_' ) then
            !this%pp = this%pp + 1
            !do while ( (c >= 'a' .and. c <= 'z') .or. c == '_' .or.  (c >= '0' .and. c <= '9') )
            !    this%pp = this%pp + 1
            !    c = this%expr(this%pp:this%pp)
            !enddo
        do j = 1, NUM_FUNCTIONS
            d = is_buildin(this, functions(j))
            if ( d > 0) exit
        enddo
        if( d > 0 ) then
            this%pp = this%pp + d  ! consume operator
            c = this%expr(this%pp:this%pp)
            if (c == '(') then
                ! check if the name matches any buildin function
                if (j > NUM_FUNCTIONS) call handle_error(this, stat%code, (/err_func_nf, pp/))

                if (DEBUG_VERBOSE) print*, 'descending into function ', this%expr(pp: this%pp-1)
                this%pp = this%pp + 1  ! consume '('
                np = 0  ! parameter count
                do
                    call l0expr(this, stat)
                    np = np + 1
                    mo(np) = stat%code(2)
                    if ( this%expr(this%pp:this%pp) /= ',' ) exit
                    this%pp = this%pp + 1  ! consume ','
                enddo

                ! check for closing parenthesis
                if (this%expr(this%pp:this%pp) /= ')') call handle_error(this, stat%code, (/err_brt_mtch, this%pp/))

                ! check call signature
                if (np /= functions(j)%n_operands) call handle_error(this, stat%code, (/err_func_ipc, pp/))

                ! matched, record instructions
                this%pp = this%pp + 1  ! consume ')'
                call write_stat( this, stat, this%pp, (/ functions(j)%id, this%rp, mo(1:np) /))
                if (DEBUG) print*, 'matched function "', this%expr(pp:this%pp-1), '" instruction = ', stat%code(1:np+2) 
            endif
        else
            call write_stat( this, stat, this%pp, (/ flag_failed /))
        endif
    end subroutine

!----------------------------------------------------------------------------------------------------------
! brt     = '(' L0EXPR ')'
!----------------------------------------------------------------------------------------------------------
    recursive subroutine bracket(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat

        ! local
        integer :: pp
        pp = this%pp
        if (this%expr(pp:pp) /= '(') then
            call write_stat( this, stat, pp, (/ flag_failed /))
            return
        endif

        this%pp = this%pp + 1  ! consume '('
        call l0expr(this, stat)
        if (stat%code(1) == flag_failed) call handle_error(this, stat%code, (/err_brt_cntt, pp+1/))

        if (this%expr(this%pp:this%pp) /= ')') call handle_error(this, stat%code, (/err_brt_mtch, this%pp/))

        call write_stat( this, stat, this%pp, (/ op_noop%id, stat%code(2) /) )

        if (DEBUG_VERBOSE) print*, 'matched parentheses "', this%expr(pp:this%pp), '"'

        this%pp = this%pp + 1  ! consume ')'

    end subroutine

!----------------------------------------------------------------------------------------------------------
! match constant defined in the constants list
!----------------------------------------------------------------------------------------------------------
    recursive subroutine constant(this, stat)
        implicit none
        type(t_parser),       intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat

        ! local
        integer :: i
        do i = 1, NUM_CONSTANT
            if ( trim(constants(i)%name) == this%expr( this%pp: this%pp + len_trim( constants(i)%name )-1 ) ) then
                this%pp = this%pp + len_trim(constants(i)%name) ! consume constant
                call write_stat( this, stat, this%pp, (/ flag_load, i + this%nvar /))
                if (DEBUG) print*, 'matched constant "', this%expr(this%pp - len_trim(constants(i)%name): this%pp-1), '"'
                exit
            endif
        enddo
    end subroutine

!----------------------------------------------------------------------------------------------------------
! setup the parser, parse input, compile and (TODO: optimize)
!----------------------------------------------------------------------------------------------------------
    subroutine setup_parser(this, str_expr, varname, ierr)
        implicit none
        type(t_parser),                         intent (inout) :: this
        character ( len = * ),               intent (in)    :: str_expr
        character ( len = * ), dimension(:), intent (in)    :: varname
        integer,                             intent (out)   :: ierr 
        
        ! local
        integer :: i, minl = huge(1), maxl = 0
        type(t_parse_result) :: stat
        ierr = -100
        call pre_processing(str_expr, this%expr, this%slen)
        this%nvar = size(varname)

        allocate( this%var_name( this%nvar ) )
        do i = 1, this%nvar
            this%var_name(i) = varname(i)
        enddo
        this%rp = 1 + NUM_CONSTANT + this%nvar
        this%pp = 1
        this%ip = 0

        allocate (this%inst_parse(ILEN, this%slen))
        allocate (this%reg_parse(this%slen/2+1 + this%rp))

        call parse(this, stat)
        ierr = stat%code(1)
        call compile_parser(this, stat)
        call cleanup_parser(this)
    end subroutine setup_parser

!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
    subroutine parse(this, stat)
        implicit none
        type(t_parser),          intent (inout) :: this
        type(t_parse_result), intent (out) :: stat
        
        !local
        if (DEBUG) print*, 'parsing ', trim(this%expr)
        call l4expr(this, stat)
        if (stat%next < len_trim(this%expr)) call handle_error(this, stat%code, (/err_eof, stat%next/))
    end subroutine parse

!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
    subroutine compile_parser(this, stat)
        implicit none
        type(t_parser),          intent (inout) :: this
        type(t_parse_result), intent (inout) :: stat

        ! local
        integer :: icnt, rcnt, i
        icnt = this%ip
        rcnt = 0
        ! remove noop/loading instructions
        if (this%rp > NUM_CONSTANT + this%nvar + 1) then
            do i = 1, this%ip
                if (this%inst_parse(1, i) == flag_load) rcnt = rcnt + 1
            enddo
        endif
        ! allocate proper size register and instruction buffer
        this%ip = this%ip - rcnt
        if (this%ip .le. 0) this%ip = 1
        this%rp = this%rp - 1
        if (DEBUG) print*, 'compiling ...  register count =', this%rp, ', instruction count =', this%ip
        allocate (this%inst(ILEN, this%ip))
        allocate (this%reg(this%rp))
        ! copy buildin constants
        do i = 1, NUM_CONSTANT 
            this%reg(this%nvar + i) = constants(i)%val
        enddo
        this%reg(NUM_CONSTANT+this%nvar+1:this%rp) = this%reg_parse(NUM_CONSTANT+this%nvar+1:this%rp)
        rcnt = 0
        this%inst(:, 1) = this%inst_parse(:, 1)
        do i = 2, icnt
            if (this%inst_parse(1, i) /= flag_load) then
                rcnt = rcnt + 1
                this%inst(:, rcnt) = this%inst_parse(:, i)
            endif
        enddo
    end subroutine compile_parser

!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
    function eval_parser(this, var)
        implicit none
        type(t_parser),     intent(inout) :: this
        real(p_k_parse), dimension(:)  :: var
        real(p_k_parse)   :: eval_parser

        !local
        integer :: i, j, k, iw, flag
        integer, dimension(ILEN) :: inst
        ! copy variable to register
        if (DEBUG) print*, '========================= INSTRUCTIONS =========================='
        this%reg(1:this%nvar) = var(1:this%nvar)
        do i = 1, this%ip
            inst = this%inst(:, i)

            if (DEBUG) print*, '>->',  inst
            if (DEBUG_RVM) then
                write(*,'(a)',advance="no") 'reg:'
                do j = 1, this%rp
                    write(*,'(a,I0,a,F0.6)',advance="no")  ' [', j, ']', this%reg(j)
                enddo
                print*, ''
                flag = 0
                do j = 1, NUM_OPRS
                    if ( inst(1) == operators(j)%id ) then
                        flag = 1
                        iw = index(operators(j)%name, ' ')
                        write(*,'(a,I0,a,F0.6,a,F0.6)',advance="no") '(*executing*)  reg[', inst(2), &
                            & '] = ', this%reg(inst(3)), operators(j)%name(1:iw), this%reg(inst(4))
                    endif
                    if ( flag == 1 ) exit
                enddo
                do j = 1, NUM_FUNCTIONS
                    if ( flag == 1 ) exit
                    if ( inst(1) == functions(j)%id ) then
                        flag = 1
                        iw = index(functions(j)%name, ' ')
                        write(*,'(a,I0,a,a,a)',advance="no") '(*executing*)  reg[', inst(2), &
                            & '] = ', functions(j)%name(1:iw), '('
                        do k = 3, functions(j)%n_operands + 1
                            write(*,'(F0.6,a)',advance="no") this%reg(inst(k)), ','
                        enddo
                        write(*,'(F0.6,a)',advance="no"), this%reg(inst(functions(j)%n_operands+2)), ')'
                    endif
                enddo
            endif

            select case(this%inst(1, i))
            case(op_plus%id)
                this%reg(inst(2)) =  this%reg(inst(3)) + this%reg(inst(4))

            case(op_minus%id)
                this%reg(inst(2)) =  this%reg(inst(3)) - this%reg(inst(4))

            case(func_neg%id)
                this%reg(inst(2)) =  - this%reg(inst(3))

            case(op_mul%id)
                this%reg(inst(2)) =  this%reg(inst(3)) * this%reg(inst(4))

            case(op_div%id)
                this%reg(inst(2)) =  this%reg(inst(3)) / this%reg(inst(4))

            case(op_mod%id)
                this%reg(inst(2)) =  mod(this%reg(inst(3)), this%reg(inst(4)))

            case(op_or%id)
                if( this%reg(inst(3)) == 1.0_p_k_parse .or. this%reg(inst(4)) == 1.0_p_k_parse ) then
                    this%reg(inst(2)) = 1.0_p_k_parse
                else
                    this%reg(inst(2)) = 0.0_p_k_parse
                endif

            case(op_and%id)
                if( this%reg(inst(3)) == 1.0_p_k_parse .and. this%reg(inst(4)) == 1.0_p_k_parse ) then
                    this%reg(inst(2)) = 1.0_p_k_parse
                else
                    this%reg(inst(2)) = 0.0_p_k_parse
                endif

            case(op_not%id)
                if( this%reg(inst(3)) == 1.0_p_k_parse ) then
                    this%reg(inst(2)) = 0.0_p_k_parse
                else
                    this%reg(inst(2)) = 1.0_p_k_parse
                endif

            case(op_le%id)
                if( this%reg(inst(3)) <= this%reg(inst(4)) ) then
                    this%reg(inst(2)) = 1.0_p_k_parse
                else
                    this%reg(inst(2)) = 0.0_p_k_parse
                endif

            case(op_lt%id)
                if( this%reg(inst(3)) < this%reg(inst(4)) ) then
                    this%reg(inst(2)) = 1.0_p_k_parse
                else
                    this%reg(inst(2)) = 0.0_p_k_parse
                endif

            case(op_ge%id)
                if( this%reg(inst(3)) >= this%reg(inst(4)) ) then
                    this%reg(inst(2)) = 1.0_p_k_parse
                else
                    this%reg(inst(2)) = 0.0_p_k_parse
                endif

            case(op_gt%id)
                if( this%reg(inst(3)) > this%reg(inst(4)) ) then
                    this%reg(inst(2)) = 1.0_p_k_parse
                else
                    this%reg(inst(2)) = 0.0_p_k_parse
                endif

            case(op_eq%id)
                if( this%reg(inst(3)) == this%reg(inst(4)) ) then
                    this%reg(inst(2)) = 1.0_p_k_parse
                else
                    this%reg(inst(2)) = 0.0_p_k_parse
                endif

            case(op_ne%id)
                if( this%reg(inst(3)) /= this%reg(inst(4)) ) then
                    this%reg(inst(2)) = 1.0_p_k_parse
                else
                    this%reg(inst(2)) = 0.0_p_k_parse
                endif

            case(op_pow%id)
                this%reg(inst(2)) = this%reg(inst(3))**int(this%reg(inst(4)))

            case(func_abs%id)
                this%reg(inst(2)) =  abs(this%reg(inst(3)))

            case(func_sinh%id)
                this%reg(inst(2)) =  sinh(this%reg(inst(3)))

            case(func_sin%id)
                this%reg(inst(2)) =  sin(this%reg(inst(3)))

            case(func_cosh%id)
                this%reg(inst(2)) =  cosh(this%reg(inst(3)))

            case(func_cos%id)
                this%reg(inst(2)) =  cos(this%reg(inst(3)))

            case(func_tanh%id)
                this%reg(inst(2)) =  tanh(this%reg(inst(3)))

            case(func_tan%id)
                this%reg(inst(2)) =  tan(this%reg(inst(3)))

            case(func_exp%id)
                this%reg(inst(2)) =  exp(this%reg(inst(3)))

            case(func_log10%id)
                this%reg(inst(2)) =  log10(this%reg(inst(3)))

            case(func_log%id)
                this%reg(inst(2)) =  log(this%reg(inst(3)))

            case(func_asin%id)
                this%reg(inst(2)) =  asin(this%reg(inst(3)))

            case(func_acos%id)
                this%reg(inst(2)) =  acos(this%reg(inst(3)))

            case(func_atan2%id)
                this%reg(inst(2)) =  atan2(this%reg(inst(3)), this%reg(inst(4)))

            case(func_atan%id)
                this%reg(inst(2)) =  atan(this%reg(inst(3)))

            case(func_sqrt%id)
                this%reg(inst(2)) =  sqrt(this%reg(inst(3)))

            case(func_if%id)
                if(this%reg(inst(3)) == 1.0) then
                    this%reg(inst(2)) = this%reg(inst(4))
                else
                    this%reg(inst(2)) = this%reg(inst(5))
                endif

            case(func_pow%id)
                this%reg(inst(2)) =  this%reg(inst(3)) ** this%reg(inst(4))

            case(func_int%id)
                this%reg(inst(2)) =  int(this%reg(inst(3)))

            case(func_nint%id)
                this%reg(inst(2)) =  nint(this%reg(inst(3)))

            case(func_ceiling%id)
                this%reg(inst(2)) =  ceiling(this%reg(inst(3)))

            case(func_floor%id)
                this%reg(inst(2)) =  floor(this%reg(inst(3)))

            case(func_modulo%id)
                this%reg(inst(2)) =  modulo(this%reg(inst(3)), this%reg(inst(4)))

            case(func_rect%id)
                if ( abs(this%reg(inst(3))) <= 0.5 ) then
                    this%reg(inst(2)) =  1.0
                else
                    this%reg(inst(2)) =  0.0
                endif

            case(func_step%id)
                if ( this%reg(inst(3)) >= 0.0 ) then
                    this%reg(inst(2)) =  1.0
                else
                    this%reg(inst(2)) =  0.0
                endif

            case(func_min3%id)
                this%reg(inst(2)) =  min(this%reg(inst(3)), this%reg(inst(4)), this%reg(inst(5)))

            case(func_min%id)
                this%reg(inst(2)) =  min(this%reg(inst(3)), this%reg(inst(4)))

            case(func_max3%id)
                this%reg(inst(2)) =  max(this%reg(inst(3)), this%reg(inst(4)), this%reg(inst(5)))

            case(func_max%id)
                this%reg(inst(2)) =  max(this%reg(inst(3)), this%reg(inst(4)))

            case default
                if (this%inst(1,i) < flag_load) call handle_error(this, inst, (/ err_runtime, 0 /))
            end select
            if ( DEBUG_RVM ) then
                write(*,'(a,F0.6)',advance="yes"), ' =', this%reg(inst(2))
            endif
        enddo
        eval_parser = this%reg(this%inst(2, i-1))
    end function

!----------------------------------------------------------------------------------------------------------
! remove all whitespaces and convert to lower case
!----------------------------------------------------------------------------------------------------------
    subroutine pre_processing( in_str, out_str, slen )
        character ( len = * ), intent(in)    :: in_str
        character ( len = * ), intent(inout) :: out_str
        integer, optional    , intent(inout) :: slen
    
        integer :: i, j
        j = 1
        do i = 1, len(in_str)
            if (in_str(i:i) /= achar(9) .and. in_str(i:i) /= ' ') then
                if ('A' <= in_str(i:i) .and. in_str(i:i) <=  'Z') then
                    out_str(j:j) = achar(iachar(in_str(i:i))+32)
                else
                    out_str(j:j) = in_str(i:i)
                endif
                j = j +1
            endif
        enddo
        if (present(slen)) slen = j
    end subroutine

!----------------------------------------------------------------------------------------------------------
! print error message
!----------------------------------------------------------------------------------------------------------
    subroutine handle_error(this, err_out, err_code, s)

        type(t_parser),     intent(in) :: this
        integer, dimension(:), intent(inout) :: err_out
        integer, dimension(2), intent(in) :: err_code
        character (len=*), optional, intent(in) :: s

        integer :: i

        select case (err_code(1))
            case ( err_opr )
                if ( present(s) ) then
                    print*, 'No operand found after ' // s
                else
                    print*, 'Missing operand'
                endif
            case ( err_func_nf )
                print*, 'Function not found'
            case ( err_func_ipc )
                print*, 'Incorrect number of function parameters'
            case ( err_brt_mtch )
                print*, 'Unbalance bracket, ")" missing'
            case ( err_brt_cntt )
                print*, 'Failed to parse contents inside parentheses'
            case ( err_var )
                print*, 'Invalid variable name'
            case ( err_num )
                print*, 'Invalid number'
            case ( err_eof )
                print*, 'Parser does not consume all characters'
            case ( err_noterm )
                print*, 'No terminal term found'
            case ( err_runtime )
                print*, 'Runtime error: unknown instruction', err_out
                stop
            case default
                print*, "Parser error"
        end select
        
        err_out(1:2) = err_code
        print*, trim(this%expr)
        write(*,*) ('-',i=1,err_code(2)-1), '^'
        print*, "error occured at position = ", err_code(2)
        stop
    end subroutine handle_error
end module

!program main
!    call test('1', (/'x1'/), (/2.1/))
!    call test('+2^2.0^3.0', (/'x1'/), (/2.1/))  ! should be 2^8 = 256
!    call test('-2 + x**y * pi', (/'x', 'y'/), (/2.0, 8./))
!    call test('1+x1+2', (/'x1'/), (/0.9/))
!    call test('if(x>0, max3(x, y, 10), atan(x-y)/10)', (/'x', 'y'/), (/1.0, 2.0/))
!    call test('x1>0', (/'x1'/), (/0.9/)) ! this is not allowed
!
!    contains
!    subroutine test( str, varname, var)
!        use m_rvm
!        character (len=*), intent(in) :: str
!        character (len=*), dimension(:) :: varname
!        real, dimension(:) :: var
!
!        !local
!        type(t_parser) :: parser
!        real(p_k_parse) :: res
!        integer :: ierr
!        call setup(parser, str, varname, ierr)
!        res = eval(parser, real(var, p_k_parse))
!        print*, 'result = ', res
!        print*, ''
!        call delete(parser)
!    end subroutine test
!
!end
!
