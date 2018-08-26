#ELSD - Ellipse and Line Segment Detector
#
#  Copyright (c) 2012 viorica patraucean (vpatrauc@gmail.com)
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program. If not, see <http://www.gnu.org/licenses/>.

DIR_F2C=../f2c/
DIR_BLAS=/opt/local/include/

elsd :	elsd.c valid_curve.c process_curve.c process_line.c write_svg.c
	cc -o elsd_io -I${DIR_F2C} -I${DIR_BLAS}  elsd_io.c valid_curve.c process_curve.c process_line.c write_svg.c ${DIR_F2C}/libf2c.a -llapack -lblas  -lm




