/********************************************************************
mkprof (MakeProfiles) - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>

#include "checkset.h"
#include "fitsarrayvv.h"

#include "main.h"





/*

Basic Strategy
==============

N: Number of threads.

QUEUE
-----

Queue of profiles that are built and ready to be added to the final
image. Each with 5 members: a pointer to the array, the position (in
the final output image of where the array should start) and the
position it should end.

The queue should be two way (with each element having a pointer to its
next and previous elements).


BUILDER THREADS
---------------

N-1 threads build profiles and fill an internal queue identical to the
global one but only for new models they have built since the last time
they added their constructs to the queue.

Upon building each profile, using pthread_mutex_trylock, they check to
see if they can lock the mutex for the global queue. If so, they add
their queue to the global one and free the mutex. Then they reset
their internal queue. If they can't get the lock, they continue adding
to their internal queue until they get the lock.


ADDER THREAD
------------

One thread is in charge of adding the separate profiles to the output
array using the global queue. Using a conditional, it waits until the
queue has some elements to start building. It starts building from the
end of the queue (note that the threads add elements from the start of
the queue). Every element it builts, it frees the array within it and
the element its self, then goes onto the next queue. If it reaches the
start of the queue and all the threads are finished, then the job is
done.


 */




















/**************************************************************/
/************           Outside function          *************/
/**************************************************************/
void
mkprof(struct mkprofparams *p)
{

}
