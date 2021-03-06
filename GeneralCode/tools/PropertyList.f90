!
! PropertyList - a Java-style Property list
!
! Property File Format is :
!
!   <Key> : <Value>  # Comment
!
! where <Key> is the name of the parameter and <Value>
! is its value.  The '#' character denotes a comment for the rest
! of the line
!
! CAVEAT: DO NOT USE TABS IN THE FILE OR READING THE DATA WILL FAIL
!
! Interface:
!
!   - LoadPropertyFile: load a property file from disk
!   - SavePropertyFile: save a property file to disk
!   - SetVal: set a property value
!   - GetVal: get a property value
!
! WARNINGS:
!   * CAVEAT: DO NOT USE TABS IN THE FILE OR READING THE DATA WILL FAIL
!   * PropertyList is not reentrant - it is only safe to use from one
!     thread of execution at a time!
!   * PropertyList has a maximum number of elements of XXX
!   * PropertyList should be implemented with something like a std::map<>
!     but Fortran lacks such niceties so I have made do with parallel
!     arrays.  If performance matters, this will need to be replaced.
!
! modification history
! --------------------
! 08apr08 bss increased PL_VAL_LEN to handle even larger tlife values.
! 24jan08 bss increased MAX_LEN to prevent truncation of long lines.
! 27nov07 bss increased PL_VAL_LEN to handle large tlife values in
!             housing application.
! 19sep07 bss fixed bug -- SetVal failed if PropertyFile is empty!
! 07sep07 bss fixed bug -- LoadPropertyFile failed to close prop file!
! 06jul07 bss written.
!

MODULE PropertyList

IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declarations

INTEGER, PARAMETER	:: PL_ERROR = -1	! ERROR return code
INTEGER, PARAMETER	:: PL_OK = 0		! Success return code

INTEGER, PARAMETER	:: PL_MAX_ELEM = 128  ! Max num properties in list
INTEGER, PARAMETER	:: PL_KEY_LEN  = 128  ! Max len of name of property
INTEGER, PARAMETER	:: PL_VAL_LEN  = 8192 ! Max len of corresponding value

CHARACTER, PARAMETER	:: PL_COMMENT_CHAR='#'	! Comment char in prop file
CHARACTER, PARAMETER	:: PL_COMMENT_SEP=':'	! Separator char in prop file

! A type to hold Property information
TYPE Property
  PRIVATE
    CHARACTER( LEN=PL_KEY_LEN ), DIMENSION( PL_MAX_ELEM ) :: keyList
    CHARACTER( LEN=PL_VAL_LEN ), DIMENSION( PL_MAX_ELEM ) :: valList
    LOGICAL		:: bInitialized = .FALSE. ! TRUE if prop data loaded
    INTEGER		:: nElem = 0		! number of elements
  
END TYPE Property

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define interface

  PUBLIC :: LoadPropertyFile
  PUBLIC :: SavePropertyFile
  PUBLIC :: SetVal
  PUBLIC :: GetVal
  PUBLIC :: isInitialized

  PRIVATE :: FindKey

CONTAINS

! LoadPropertyFile - loads property information from a file
!
! Params:
!   prop_	 - the Property object
!   szFilename_  - file to load
!   retval - PL_ERROR on failure, PL_OK otherwise

  INTEGER FUNCTION LoadPropertyFile( prop_, szFilename_ )
    TYPE( Property ) :: prop_
    CHARACTER(LEN=*) :: szFilename_ 

    INTEGER, PARAMETER		:: MAX_LEN = PL_KEY_LEN + PL_VAL_LEN + 1	! maximum length of line in file
    INTEGER, PARAMETER		:: nInUnit = 11

    INTEGER			:: ios	! I/O Status
    CHARACTER(LEN=MAX_LEN)		:: szLine
    CHARACTER(LEN=PL_KEY_LEN)	:: szKey
    CHARACTER(LEN=PL_VAL_LEN)	:: szVal
    INTEGER			:: ixPosCom	! comment position in string
    INTEGER			:: ixPosSep   ! key/val sep position in string
    INTEGER			:: ix		! loop index

    OPEN(UNIT=nInUnit, FILE=szFilename_, STATUS='OLD', &
         ACCESS='SEQUENTIAL', ACTION='READ', IOSTAT=ios)

    IF( ios /= 0 ) THEN
      LoadPropertyFile = PL_ERROR
      RETURN
    END IF

    ix = 0 
    DO 
      READ(UNIT=nInUnit, FMT='(A)', END=1000 ) szLine

      IF( ios /= 0 ) THEN
        LoadPropertyFile = PL_ERROR
        EXIT
      ENDIF

      ixPosCom = INDEX( szLine, PL_COMMENT_CHAR )
      ixPosSep = INDEX( szLine, PL_COMMENT_SEP )
      
      ! skip lines which are comments or lack a separator
      IF( (ixPosCom == 1) .OR.  (ixPosSep == 0) ) THEN
        CYCLE
      END IF 

      ix = ix + 1

      ! Calculate end of line
      IF( ixPosCom == 0 ) THEN
        ixPosCom = MAX_LEN
      ELSE
        ixPosCom = ixPosCom - 1
      ENDIF

      szKey = szLine( 1:ixPosSep-1 )
      szVal = szLine( ixPosSep+1 : ixPosCom )

      szKey = TRIM( ADJUSTL( szKey ) )
      prop_%keyList(ix) = szKey
      prop_%keyList(ix) = TRIM( ADJUSTL( szKey ) )
      prop_%valList(ix) = TRIM( ADJUSTL( szVal ) )

    END DO

1000 LoadPropertyFile = PL_OK
    prop_%bInitialized = .TRUE.     
    prop_%nElem = ix

    CLOSE( nInUnit )

  END FUNCTION

! SavePropertyFile - saves property information to a file
!
! Params:
!   prop_	 - the Property object
!   szFilename_  - file to load
!   retval - PL_ERROR on failure, PL_OK otherwise

  INTEGER FUNCTION SavePropertyFile( prop_, szFilename_ )
    TYPE( Property ) :: prop_
    CHARACTER(LEN=*) :: szFilename_ 
    
    INTEGER, PARAMETER		:: MAX_LEN = PL_KEY_LEN + PL_VAL_LEN + 1	! maximum length of line in file
    INTEGER, PARAMETER		:: nOutUnit = 12

    INTEGER			:: ios	! I/O Status
    INTEGER			:: ix		! loop index
    
    SavePropertyFile = PL_ERROR

    IF( .NOT. prop_%bInitialized ) THEN
      RETURN
    END IF

    OPEN(UNIT=nOutUnit, FILE=szFilename_, STATUS='REPLACE', &
         ACCESS='SEQUENTIAL', ACTION='WRITE', IOSTAT=ios) 

    IF( ios /= 0 ) THEN
      SavePropertyFile = PL_ERROR
      RETURN
    END IF

    DO ix = 1, prop_%nElem
      WRITE( nOutUnit, '(3A)' ) TRIM( prop_%keyList(ix) ), ' : ', &
             TRIM( prop_%valList(ix) )
    END DO

    CLOSE( UNIT=nOutUnit )

    SavePropertyFile = PL_OK
  END FUNCTION

! SetVal - changes the value associated with a give key
!
! If szKey is not found, a new value is created
!
! Params:
!   prop_	 - the Property object
!   szKey_ - the key whose value should be modified
!   szVal_ - the value to store
!   retval - PL_ERROR/PL_OK

  INTEGER FUNCTION SetVal( prop_, szKey_, szVal_ )
    TYPE( Property ) 			:: prop_
    CHARACTER(LEN=*), INTENT(IN)	:: szKey_
    CHARACTER(LEN=*), INTENT(IN) 	:: szVal_ 
    INTEGER				:: ixKey	! index of szKey_

    SetVal = PL_ERROR

    IF( .NOT. prop_%bInitialized ) THEN
      prop_%bInitialized = .TRUE.
!      SetVal = PL_ERROR
!      RETURN
    ENDIF

    ixKey = FindKey( prop_, szKey_ )
    IF( ixKey == PL_ERROR ) THEN
      prop_%nElem = prop_%nElem + 1
      ixKey = prop_%nElem
      prop_%keyList(ixKey) = szKey_
    ENDIF

    prop_%valList(ixKey) = szVal_ 
    SetVal = PL_OK
    
  END FUNCTION SetVal

! GetVal - gets the value associated with a give key
!
! Params:
!   prop_	 - the Property object
!   szKey_ - the key whose value is desired
!   szVal_ - the value of the key
!   retval - PL_ERROR/PL_OK

  INTEGER FUNCTION GetVal( prop_, szKey_, szVal_ )
    TYPE( Property ) 			:: prop_
    CHARACTER(LEN=*), INTENT(IN) 	:: szKey_
    CHARACTER(LEN=*), INTENT(OUT) 	:: szVal_ 
    INTEGER				:: ixKey	! index of szKey_

    GetVal = PL_ERROR
    ixKey = FindKey( prop_, szKey_ )
    IF( ixKey /= PL_ERROR ) THEN
      szVal_ = prop_%valList(ixKey)
      GetVal = PL_OK
    ENDIF
    
  END FUNCTION GetVal

! isInitialized - returns .TRUE. iff Property has been initialized
!
! Params:
!   prop_ - Property list to check
!   retVal - .TRUE./.FALSE.

  FUNCTION isInitialized( prop_ )
    LOGICAL :: isInitialized
    TYPE( Property ) :: prop_

    isInitialized = prop_%bInitialized

  END FUNCTION isInitialized

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Private Procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! FindKey - find the index where a key is stored
!
! Params
!   prop_	 - the Property object
!   szKey_ - the key to find
!   retval - the index of szKey_ or PL_ERROR

  INTEGER FUNCTION FindKey( prop_, szKey_ )
    TYPE( Property ) 			:: prop_
    CHARACTER(LEN=*), INTENT(IN) 	:: szKey_

    INTEGER 		:: ix		! loop index
    
    FindKey = PL_ERROR			! Assume failure

    DO ix = 1, prop_%nElem
      ! exit loop when key is found
      IF( szKey_ == prop_%keyList(ix) ) THEN
        FindKey = ix
        EXIT
      END IF
    END DO
  END FUNCTION FindKey

END MODULE PropertyList
