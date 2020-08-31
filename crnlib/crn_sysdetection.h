#pragma once

// OS Detection

/*
   The operating system, must be one of: (CRN_OS_x)
     DARWIN   - Any Darwin system (macOS, iOS, watchOS, tvOS)
     MACOS    - macOS
     IOS      - iOS
     WATCHOS  - watchOS
     TVOS     - tvOS
     WIN32    - Win32 (Windows 2000/XP/Vista/7 and Windows Server 2003/2008)
     CYGWIN   - Cygwin
     SOLARIS  - Sun Solaris
     LINUX    - Linux [has variants]
     FREEBSD  - FreeBSD [has variants]
     NETBSD   - NetBSD
     OPENBSD  - OpenBSD
     BSD4     - Any BSD 4.4 system
     HURD     - GNU Hurd
     UNIX     - Any UNIX BSD/SYSV system
     ANDROID  - Android platform
     WEBOS    - LG WebOS
   The following operating systems have variants:
     LINUX    - both CRN_OS_LINUX and CRN_OS_ANDROID are defined when building for Android
              - only CRN_OS_LINUX is defined if building for other Linux systems
     MACOS    - both CRN_OS_BSD4 and CRN_OS_IOS are defined when building for iOS
              - both CRN_OS_BSD4 and CRN_OS_MACOS are defined when building for macOS
     FREEBSD  - CRN_OS_FREEBSD is defined only when building for FreeBSD with a BSD userland
              - CRN_OS_FREEBSD_KERNEL is always defined on FreeBSD, even if the userland is from GNU
*/

#if defined(__APPLE__) && (defined(__GNUC__) || defined(__xlC__) || defined(__xlc__))
#  include <TargetConditionals.h>
#  if defined(TARGET_OS_MAC) && TARGET_OS_MAC
#    define CRN_OS_DARWIN
#    define CRN_OS_BSD4
#    ifdef __LP64__
#      define CRN_OS_DARWIN64
#    else
#      define CRN_OS_DARWIN32
#    endif
#    if defined(TARGET_OS_IPHONE) && TARGET_OS_IPHONE
#      if defined(TARGET_OS_WATCH) && TARGET_OS_WATCH
#        define CRN_OS_WATCHOS
#      elif defined(TARGET_OS_TV) && TARGET_OS_TV
#        define CRN_OS_TVOS
#      else
#        // TARGET_OS_IOS is only available in newer SDKs,
#        // so assume any other iOS-based platform is iOS for now
#        define CRN_OS_IOS
#      endif
#    else
#      // TARGET_OS_OSX is only available in newer SDKs,
#      // so assume any non iOS-based platform is macOS for now
#      define CRN_OS_MACOS
#    endif
#  else
#    error "crn doesn't support this Apple platform"
#  endif
#elif defined(__WEBOS__)
#  define CRN_OS_WEBOS
#elif defined(__ANDROID__) || defined(ANDROID)
#  define CRN_OS_ANDROID
#elif defined(__CYGWIN__)
#  define CRN_OS_CYGWIN
#elif !defined(SAG_COM) && (!defined(WINAPI_FAMILY) || WINAPI_FAMILY==WINAPI_FAMILY_DESKTOP_APP) && (defined(WIN64) || defined(_WIN64) || defined(__WIN64__))
#  define CRN_OS_WIN32
#  define CRN_OS_WIN64
#elif !defined(SAG_COM) && (defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__))
#  define CRN_OS_WIN32
#elif defined(__sun) || defined(sun)
#  define CRN_OS_SOLARIS
#elif defined(__native_client__)
#  define CRN_OS_NACL
#elif defined(__EMSCRIPTEN__)
#  define CRN_OS_WASM
#elif defined(__linux__) || defined(__linux)
#  define CRN_OS_LINUX
#elif defined(__FreeBSD__) || defined(__DragonFly__) || defined(__FreeBSD_kernel__)
#  ifndef __FreeBSD_kernel__
#    define CRN_OS_FREEBSD
#  endif
#  define CRN_OS_FREEBSD_KERNEL
#  define CRN_OS_BSD4
#elif defined(__NetBSD__)
#  define CRN_OS_NETBSD
#  define CRN_OS_BSD4
#elif defined(__OpenBSD__)
#  define CRN_OS_OPENBSD
#  define CRN_OS_BSD4
#elif defined(__GNU__)
#  define CRN_OS_HURD
#else
#  error "crn doesn't support this OS"
#endif

#if defined(CRN_OS_WIN32) || defined(CRN_OS_WIN64)
#  define CRN_OS_WIN
#endif

#if defined(CRN_OS_WIN)
#  undef CRN_OS_UNIX
#elif !defined(CRN_OS_UNIX)
#  define CRN_OS_UNIX
#endif

#ifdef CRN_OS_DARWIN
#define CRN_OS_MAC
#endif
#ifdef CRN_OS_DARWIN32
#define CRN_OS_MAC32
#endif
#ifdef CRN_OS_DARWIN64
#define CRN_OS_MAC64
#endif
#ifdef CRN_OS_MACOS
#define CRN_OS_MACX
#define CRN_OS_OSX
#endif

// Compiler detection



