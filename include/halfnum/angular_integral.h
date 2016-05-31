/***********
This file is part of libFireDeamon.

Copyright (C) 2015,2016 by Torsten Sachse

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
#ifndef HALFNUM_ANGULAR_INTEGRALS_H
#define HALFNUM_ANGULAR_INTEGRALS_H

#define LMAXP1 6
class AngInt{
    private:
        double* m_integrals;
    public:
        AngInt();
        double GetInt(unsigned int lambda, int mu, unsigned int i, unsigned int j, unsigned int k) const;
        ~AngInt();
};
#endif //HALFNUM_ANGULAR_INTEGRALS_H
