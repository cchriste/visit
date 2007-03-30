#ifndef DBATTS_EXPORTS_H
#define DBATTS_EXPORTS_H

#if defined(_WIN32)
#ifdef DBATTS_EXPORTS
#define DBATTS_API __declspec(dllexport)
#else
#define DBATTS_API __declspec(dllimport)
#endif
#ifdef USING_MSVC6
// Turn off warning about lack of DLL interface
#pragma warning(disable:4251)
// Turn off warning non-dll class is base for dll-interface class.
#pragma warning(disable:4275)
// Turn off warning about identifier truncation
#pragma warning(disable:4786)
#endif
#else
#define DBATTS_API
#endif

#endif
