#pragma once

#include "../Core/DCoordinates3.h"

namespace cagd
{
    namespace spiral_on_cone
    {
        extern GLdouble u_min, u_max;

        DCoordinate3 d0(GLdouble);
        DCoordinate3 d1(GLdouble);
        DCoordinate3 d2(GLdouble);
    }
    namespace torus {
		extern GLdouble u_min, u_max;

		DCoordinate3 d0(GLdouble);
		DCoordinate3 d1(GLdouble);
		DCoordinate3 d2(GLdouble);
    }
	namespace lissajous {
		extern GLdouble u_min, u_max;

		DCoordinate3 d0(GLdouble);
		DCoordinate3 d1(GLdouble);
		DCoordinate3 d2(GLdouble);
	}
	namespace hypo {
		extern GLdouble u_min, u_max;

		DCoordinate3 d0(GLdouble);
		DCoordinate3 d1(GLdouble);
		DCoordinate3 d2(GLdouble);
	}
	namespace cyclo {
		extern GLdouble u_min, u_max;

		DCoordinate3 d0(GLdouble);
		DCoordinate3 d1(GLdouble);
		DCoordinate3 d2(GLdouble);
	}
	namespace ellipse {
		extern GLdouble u_min, u_max;

		DCoordinate3 d0(GLdouble);
        DCoordinate3 d1(GLdouble);
        DCoordinate3 d2(GLdouble);
	}
	namespace rose {
		extern GLdouble u_min, u_max;

		DCoordinate3 d0(GLdouble);
		DCoordinate3 d1(GLdouble);
		DCoordinate3 d2(GLdouble);
	}

	//Surface

	namespace first{
		extern GLdouble u_min, u_max, v_min, v_max;

		DCoordinate3 d00(GLdouble, GLdouble);
		DCoordinate3 d10(GLdouble, GLdouble);
		DCoordinate3 d01(GLdouble, GLdouble);
	}

	namespace rose_surface {
		extern GLdouble u_min, u_max, v_min, v_max;

		DCoordinate3 d00(GLdouble, GLdouble);
		DCoordinate3 d10(GLdouble, GLdouble);
		DCoordinate3 d01(GLdouble, GLdouble);
	}

	namespace third {
		extern GLdouble u_min, u_max, v_min, v_max;

		DCoordinate3 d00(GLdouble, GLdouble);
		DCoordinate3 d10(GLdouble, GLdouble);
		DCoordinate3 d01(GLdouble, GLdouble);
	}

	namespace fourth {
		extern GLdouble u_min, u_max, v_min, v_max;
		DCoordinate3 d00(GLdouble, GLdouble);
		DCoordinate3 d10(GLdouble, GLdouble);
		DCoordinate3 d01(GLdouble, GLdouble);

	}

	namespace fifth {
		extern GLdouble u_min, u_max, v_min, v_max;
		DCoordinate3 d00(GLdouble, GLdouble);
		DCoordinate3 d10(GLdouble, GLdouble);
		DCoordinate3 d01(GLdouble, GLdouble);
	}

	namespace sixth {
		extern GLdouble u_min, u_max, v_min, v_max;
		DCoordinate3 d00(GLdouble, GLdouble);
		DCoordinate3 d10(GLdouble, GLdouble);
		DCoordinate3 d01(GLdouble, GLdouble);
	}
    namespace dini {
        extern GLdouble u_min, u_max, v_min, v_max;
        DCoordinate3 d00(GLdouble, GLdouble);
        DCoordinate3 d10(GLdouble, GLdouble);
        DCoordinate3 d01(GLdouble, GLdouble);
    }
    namespace torus_surf {
        extern GLdouble u_min, u_max, v_min, v_max;
        DCoordinate3 d00(GLdouble, GLdouble);
        DCoordinate3 d10(GLdouble, GLdouble);
        DCoordinate3 d01(GLdouble, GLdouble);
    }
}
