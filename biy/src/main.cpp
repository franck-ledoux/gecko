// biy ("Block It Yourself"): loads a .msh file as a gecko GeomModel, shows it in an interactive 3D
// window, and lets you build a Blocking against it — from the on-screen panel, from mouse drags on
// the block's corners, or from the Python console on stdin. See docs/user-guide/biy.md.
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

#include <polyscope/options.h>
#include <polyscope/polyscope.h>

#include "BiyApp.h"
#include <gecko/app/BlockingFacade.h>
#include "PythonConsole.h"

int main(int argc, char **argv) {
    using gecko::app::BlockingFacade;

    if (argc < 2 || argc > 3) {
        std::cerr << "usage: biy <model.msh> [order]\n"
                  << "  order  block edge order the session starts at, " << BlockingFacade::MIN_DEGREE
                  << " (straight) or higher; default 3. Changeable from the panel.\n";
        return 1;
    }

    int order = 3;
    if (argc == 3) {
        try {
            std::size_t consumed = 0;
            order = std::stoi(argv[2], &consumed);
            if (consumed != std::string(argv[2]).size()) throw std::invalid_argument("trailing characters");
        } catch (const std::exception &) {
            std::cerr << "biy: order must be a whole number, got '" << argv[2] << "'\n";
            return 1;
        }
        if (order < BlockingFacade::MIN_DEGREE || order > gecko::biy::BiyApp::MAX_ORDER) {
            std::cerr << "biy: order must be between " << BlockingFacade::MIN_DEGREE << " and "
                      << gecko::biy::BiyApp::MAX_ORDER << ", got " << order << "\n";
            return 1;
        }
    }

    try {
        // Set before init(): the GLFW window is created with this as its title, so it is too late
        // to change once Polyscope is up.
        polyscope::options::programName = "BIY";
        polyscope::init();

        gecko::biy::BiyApp app(argv[1], order);
        polyscope::state::userCallback = [&app] { app.per_frame(); };

        gecko::biy::PythonConsole console(app);

        // Verification aid: render a few frames, save what the window shows, and exit — lets the
        // rendering path be checked without a human at the keyboard.
        if (const char *shot = std::getenv("BIY_SCREENSHOT")) {
            for (int i = 0; i < 5; ++i)
                polyscope::frameTick();
            polyscope::screenshot(shot);
            console.stop();
            std::cout << "Wrote screenshot to " << shot << " — press Enter to exit.\n" << std::flush;
            return 0;
        }

        // frameTick() instead of show(): it renders one frame and returns, leaving this loop free to
        // also watch for the console asking to quit.
        while (!polyscope::windowRequestsClose() && console.running()) {
            polyscope::frameTick();
        }

        if (console.running()) {
            std::cout << "\nWindow closed — press Enter to exit.\n" << std::flush;
            console.stop();
        }
    } catch (const std::exception &e) {
        std::cerr << "biy: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
