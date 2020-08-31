/****************************************************************************
**
** Copyright (C) 2014 Digia Plc and/or its subsidiary(-ies).
** Copyright (C) 2012 Intel Corporation
** Contact: http://www.qt-project.org/legal
**
** This file is part of the QtCore module of the Qt Toolkit.
**
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and Digia. For licensing terms and
** conditions see http://qt.digia.com/licensing. For further information
** use the contact form at http://qt.digia.com/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 2.1 or version 3 as published by the Free
** Software Foundation and appearing in the file LICENSE.LGPLv21 and
** LICENSE.LGPLv3 included in the packaging of this file. Please review the
** following information to ensure the GNU Lesser General Public License
** requirements will be met: https://www.gnu.org/licenses/lgpl.html and
** http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** In addition, as a special exception, Digia gives you certain additional
** rights. These rights are described in the Digia Qt LGPL Exception
** version 1.1, included in the file LGPL_EXCEPTION.txt in this package.
**
****************************************************************************/

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

/*
   The compiler, must be one of: (CRN_CC_x)
     SYM      - Digital Mars C/C++ (used to be Symantec C++)
     MSVC     - Microsoft Visual C/C++, Intel C++ for Windows
     BOR      - Borland/Turbo C++
     WAT      - Watcom C++
     GNU      - GNU C++
     COMEAU   - Comeau C++
     EDG      - Edison Design Group C++
     OC       - CenterLine C++
     SUN      - Forte Developer, or Sun Studio C++
     MIPS     - MIPSpro C++
     DEC      - DEC C++
     HPACC    - HP aC++
     USLC     - SCO OUDK and UDK
     CDS      - Reliant C++
     KAI      - KAI C++
     INTEL    - Intel C++ for Linux, Intel C++ for Windows
     HIGHC    - MetaWare High C/C++
     PGI      - Portland Group C++
     GHS      - Green Hills Optimizing C++ Compilers
     RVCT     - ARM Realview Compiler Suite
     CLANG    - C++ front-end for the LLVM compiler
   Should be sorted most to least authoritative.
*/

#if defined(__DMC__) || defined(__SC__)
#  define CRN_CC_SYM
#elif defined(_MSC_VER)
#  define CRN_CC_MSVC
#  define CRN_CC_MSVC_VER (_MSC_VER)
#  define CRN_CC_MSVC_NET
#  if defined(__INTEL_COMPILER)
#    define CRN_CC_INTEL
#  endif
#elif defined(__BORLANDC__) || defined(__TURBOC__)
#  define CRN_CC_BOR
#elif defined(__WATCOMC__)
#  define CRN_CC_WAT
#elif defined(__ARMCC__) || defined(__CC_ARM)
#  define CRN_CC_RVCT
#elif defined(__GNUC__)
#  define CRN_CC_GNU
#  define CRN_CC_GNU_VER (__GNUC__ * 100 + __GNUC_MINOR__)
#  if defined(__MINGW32__)
#    define CRN_CC_MINGW
#  endif
#  if defined(__INTEL_COMPILER) /* Intel C++ also masquerades as GCC */
#    define CRN_CC_INTEL (__INTEL_COMPILER)
#    ifdef __clang__ /* Intel C++ masquerades as Clang masquerading as GCC */
#      define CRN_CC_CLANG
#      define CRN_CC_CLANG_VER    305
#    endif
#  elif defined(__clang__) /* Clang also masquerades as GCC */
#    define CRN_CC_CLANG
#    if defined(__apple_build_version__) /* http://en.wikipedia.org/wiki/Xcode#Toolchain_versions */
#      if __apple_build_version__ >= 8020041
#        define CRN_CC_CLANG_VER 309
#      elif __apple_build_version__ >= 8000038
#        define CRN_CC_CLANG_VER 308
#      elif __apple_build_version__ >= 7000053
#        define CRN_CC_CLANG_VER 306
#      elif __apple_build_version__ >= 6000051
#        define CRN_CC_CLANG_VER 305
#      elif __apple_build_version__ >= 5030038
#        define CRN_CC_CLANG_VER 304
#      elif __apple_build_version__ >= 5000275
#        define CRN_CC_CLANG_VER 303
#      elif __apple_build_version__ >= 4250024
#        define CRN_CC_CLANG_VER 302
#      elif __apple_build_version__ >= 3180045
#        define CRN_CC_CLANG_VER 301
#      elif __apple_build_version__ >= 2111001
#        define CRN_CC_CLANG_VER 300
#      else
#        warning "Unknown Apple Clang version"
#      endif
#    else
#      define CRN_CC_CLANG_VER ((__clang_major__ * 100) + __clang_minor__)
#    endif
#  endif
#elif defined(__xlC__)
#  define CRN_CC_XLC
#elif defined(__DECCXX) || defined(__DECC)
#  define CRN_CC_DEC
#  if defined(__EDG__)
#    define CRN_CC_EDG
#  endif
#elif defined(__PGI)
#  define CRN_CC_PGI
#  if defined(__EDG__)
#    define CRN_CC_EDG
#  endif
#elif !defined(CRN_OS_HPUX) && (defined(__EDG) || defined(__EDG__))
#  define CRN_CC_EDG
#  if defined(__COMO__)
#    define CRN_CC_COMEAU
#  elif defined(__KCC)
#    define CRN_CC_KAI
#  elif defined(__INTEL_COMPILER)
#    define CRN_CC_INTEL (__INTEL_COMPILER)
#  elif defined(__ghs)
#    define CRN_CC_GHS
#  elif defined(__DCC__)
#    define CRN_CC_DIAB
#  elif defined(__USLC__) && defined(__SCO_VERSION__)
#    define CRN_CC_USLC
#  elif defined(CENTERLINE_CLPP) || defined(OBJECTCENTER)
#    define CRN_CC_OC
#  elif defined(sinix)
#    define CRN_CC_CDS
#  elif defined(__sgi)
#    define CRN_CC_MIPS
#  endif
#elif defined(_DIAB_TOOL)
#  define CRN_CC_DIAB
#elif defined(__HIGHC__)
#  define CRN_CC_HIGHC
#elif defined(__SUNPRO_CC) || defined(__SUNPRO_C)
#  define CRN_CC_SUN
#elif defined(sinix)
#  define CRN_CC_EDG
#  define CRN_CC_CDS
#elif defined(CRN_OS_HPUX) && defined(__HP_aCC)
#  define CRN_CC_HPACC
#else
#  warning "Cannot detect compiler"
#endif
