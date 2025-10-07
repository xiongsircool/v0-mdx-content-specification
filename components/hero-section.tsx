"use client"

import { motion } from 'framer-motion'
import { useEffect, useRef } from 'react'
import { gsap } from 'gsap'
import {
  Users,
  Target,
  BookOpen,
  ChevronDown,
  Code
} from 'lucide-react'
import { GlitchText, WaveText, FloatingDna, ParticleBackground } from './animations'
import Link from 'next/link'

export function HeroSection() {
  const heroRef = useRef<HTMLDivElement>(null)

  useEffect(() => {
    if (!heroRef.current) return

    const tl = gsap.timeline()

    tl.fromTo('.hero-title',
      { opacity: 0, y: 100 },
      { opacity: 1, y: 0, duration: 1, ease: 'power3.out' }
    )
    .fromTo('.hero-subtitle',
      { opacity: 0, y: 50 },
      { opacity: 1, y: 0, duration: 0.8, ease: 'power3.out' },
      '-=0.5'
    )
    .fromTo('.hero-description',
      { opacity: 0, y: 30 },
      { opacity: 1, y: 0, duration: 0.6, ease: 'power3.out' },
      '-=0.3'
    )
    .fromTo('.hero-buttons',
      { opacity: 0, y: 30 },
      { opacity: 1, y: 0, duration: 0.6, ease: 'power3.out' },
      '-=0.2'
    )

  }, [])

  return (
    <section
      ref={heroRef}
      className="relative min-h-screen flex items-center justify-center overflow-hidden bg-gradient-to-br from-gray-900 via-blue-900 to-purple-900"
    >
      <ParticleBackground className="opacity-30" />
      <FloatingDna className="opacity-20" />

  
      {/* Hero content */}
      <div className="relative z-20 max-w-6xl mx-auto px-4 text-center">
        <div className="mb-8">
          <motion.div
            className="inline-flex items-center justify-center w-24 h-24 bg-gradient-to-br from-blue-600 via-purple-600 to-pink-600 rounded-full shadow-2xl mb-8"
            animate={{
              rotate: [0, 360],
            }}
            transition={{
              duration: 20,
              repeat: Infinity,
              ease: "linear",
            }}
          >
            <Target className="w-12 h-12 text-white" />
          </motion.div>
        </div>

        <h1 className="hero-title text-6xl md:text-8xl lg:text-9xl font-bold mb-6">
          <GlitchText
            text="SBC"
            className="text-transparent bg-clip-text bg-gradient-to-r from-blue-400 via-purple-400 to-pink-400"
          />
        </h1>

        <div className="hero-subtitle mb-8">
          <WaveText
            text="上海生物信息学中心"
            className="text-2xl md:text-4xl font-bold text-white"
          />
        </div>

        <div className="hero-subtitle mb-8">
          <WaveText
            text="生物信息爱好者的聚集地"
            className="text-xl md:text-2xl text-blue-200"
          />
        </div>

        <motion.p
          className="hero-description text-lg md:text-xl text-blue-100 max-w-3xl mx-auto mb-12 leading-relaxed"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ duration: 1, delay: 1.5 }}
        >
          探索生命科学与代码的奇妙世界，<br />
          与志同道合的伙伴一起创造未来
        </motion.p>

        <div className="hero-buttons flex flex-col sm:flex-row gap-6 justify-center items-center mb-16">
          <motion.div
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
          >
            <Link
              href="/announcements"
              className="group relative px-12 py-6 bg-gradient-to-r from-blue-600 to-purple-600 text-white font-bold text-lg rounded-full shadow-2xl hover:shadow-blue-500/25 transition-all duration-300 inline-flex items-center gap-3"
            >
              <span className="relative z-10 flex items-center gap-3">
                <Users className="w-6 h-6" />
                加入我们
              </span>
              <div className="absolute inset-0 bg-gradient-to-r from-blue-600 to-purple-600 rounded-full opacity-0 group-hover:opacity-100 transition-opacity duration-300 blur-lg" />
            </Link>
          </motion.div>

          <motion.div
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
          >
            <Link
              href="/about"
              className="group relative px-12 py-6 bg-transparent border-2 border-blue-400 text-blue-300 font-bold text-lg rounded-full hover:bg-blue-400/10 transition-all duration-300 inline-flex items-center gap-3"
            >
              <span className="relative z-10 flex items-center gap-3">
                <BookOpen className="w-6 h-6" />
                了解我们
              </span>
              <div className="absolute inset-0 bg-blue-400/20 rounded-full opacity-0 group-hover:opacity-100 transition-opacity duration-300 blur-lg" />
            </Link>
          </motion.div>
        </div>

        {/* Scroll indicator */}
        <motion.div
          className="absolute bottom-8 left-1/2 transform -translate-x-1/2"
          animate={{ y: [0, 10, 0] }}
          transition={{ duration: 2, repeat: Infinity, ease: "easeInOut" }}
        >
          <ChevronDown className="w-8 h-8 text-blue-400" />
        </motion.div>
      </div>

      {/* Background gradient orbs */}
      <div className="absolute top-0 left-0 w-96 h-96 bg-blue-600 rounded-full mix-blend-multiply filter blur-3xl opacity-20 animate-pulse" />
      <div className="absolute top-0 right-0 w-96 h-96 bg-purple-600 rounded-full mix-blend-multiply filter blur-3xl opacity-20 animate-pulse [animation-delay:2s]" />
      <div className="absolute bottom-0 left-1/2 w-96 h-96 bg-pink-600 rounded-full mix-blend-multiply filter blur-3xl opacity-20 animate-pulse [animation-delay:4s]" />
    </section>
  )
}