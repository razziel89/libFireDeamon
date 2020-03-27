/***********
This file is part of libFireDeamon.

Copyright (C) 2016 by Torsten Sachse

libFireDeamon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libFireDeamon is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with libFireDeamon.  If not, see <http://www.gnu.org/licenses/>.
***********/
#include <FireDeamon/core/set_procname.h>
#include <string.h>
#include <string>
#define NAMEMAXBYTES 16
// PR_SET_NAME
#include <sys/prctl.h>
void set_procname(std::string newname) {
  char *c = (char *)malloc(NAMEMAXBYTES * sizeof(char));
  memset(c, '\0', NAMEMAXBYTES * sizeof(char));
  // Without the -1, the string is not guaranteed to end in a NULL character
  strncpy(c, newname.c_str(), (NAMEMAXBYTES - 1) * sizeof(char));
  prctl(PR_SET_NAME, c, 0, 0, 0);
}
