        subroutine datime(day)
   
        character*24 day

#ifdef PC
cDEC$ IF DEFINED (PC)
        character*24 fdate
        external fdate

        day=fdate()
cDEC$ ELSE
#elif (AIX)
        external fdate_
        call fdate_(day)
#else
        character*24 fdate

        day=fdate()
#endif
CDEC$ END IF

        end

        subroutine ddtime(tarray)
  
#ifdef PC
cDEC$ IF DEFINED (PC)
cDEC$ ELSE
#elif (AIX)
#define dtime dtime_
#endif
CDEC$ END IF

        real*4 tarray(2), sum

        call dtime(tarray, sum)

        end

        subroutine eetime(tarray)

#ifdef PC
cDEC$ IF DEFINED (PC)
cDEC$ ELSE
#elif (AIX)
#define etime etime_
#endif
CDEC$ END IF
 
        real*4 tarray(2), sum

        call etime(tarray, sum)

        end

