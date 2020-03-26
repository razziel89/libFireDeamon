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
/**
 * \file
 * \brief Function to set the name of the current process
 */
#ifndef SET_PROCNAME_H
#define SET_PROCNAME_H
// PR_SET_NAME
//#include <sys/prctl.h>
#include <string>
/**
 * \brief Set the process name
 *
 * \param newname std::string - the new name to be used
 * \param argv    char**      - the argunent vector array
 * \return whether or not setting the name succeeded
 */
void set_procname(std::string newname);

#endif // SET_PROCNAME_H
