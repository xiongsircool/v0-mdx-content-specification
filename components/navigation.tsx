"use client"

import { useState } from "react"
import Link from "next/link"
import { Button } from "@/components/ui/button"
import { Sheet, SheetContent, SheetTrigger, SheetTitle } from "@/components/ui/sheet"
import { Menu, X, BookOpen, Megaphone, GraduationCap, Calendar, Target } from "lucide-react"
import { SimpleGlobalSearch } from "@/components/simple-global-search"
import Image from "next/image"

const navigation = [
  { name: "首页", href: "/" },
  { name: "活动公告", href: "/announcements", icon: Megaphone },
  { name: "技术分享", href: "/posts", icon: BookOpen },
  { name: "学习资料", href: "/learn", icon: GraduationCap },
  { name: "关于我们", href: "/about", icon: Target },
  { name: "会议", href: "/meetings", icon: Calendar },
]

export function Navigation() {
  const [isOpen, setIsOpen] = useState(false)

  return (
    <header className="sticky top-0 z-50 w-full border-b bg-background/95 backdrop-blur supports-[backdrop-filter]:bg-background/60">
      <div className="container mx-auto px-4">
        <div className="flex h-16 items-center justify-between">
          {/* Logo */}
          <Link href="/" className="flex items-center gap-3">
            <div className="relative h-10 w-10" style={{ position: 'relative' }}>
              <Image
                src="/logo-nav.png"
                alt="SBC Logo"
                fill
                className="object-contain"
                priority
              />
            </div>
            <div className="hidden sm:block">
              <h1 className="text-lg font-bold">SBC</h1>
              <p className="text-xs text-muted-foreground">生信爱好者聚集地</p>
            </div>
          </Link>

          {/* Desktop Navigation */}
          <div className="hidden md:flex items-center gap-6">
            <nav className="flex items-center gap-6">
              {navigation.map((item) => {
                const Icon = item.icon
                return (
                  <Link
                    key={item.name}
                    href={item.href}
                    className="flex items-center gap-2 text-sm font-medium text-muted-foreground transition-colors hover:text-foreground"
                  >
                    {Icon && <Icon className="h-4 w-4" />}
                    {item.name}
                  </Link>
                )
              })}
            </nav>
            <SimpleGlobalSearch />
          </div>

          {/* Mobile Navigation */}
          <Sheet open={isOpen} onOpenChange={setIsOpen}>
            <SheetTrigger asChild className="md:hidden">
              <Button variant="ghost" size="icon">
                <Menu className="h-5 w-5" />
                <span className="sr-only">打开菜单</span>
              </Button>
            </SheetTrigger>
            <SheetContent side="right" className="w-80">
              <SheetTitle className="sr-only">导航菜单</SheetTitle>
              <div className="flex items-center justify-between mb-8">
                <div className="flex items-center gap-3">
                  <div className="relative h-8 w-8" style={{ position: 'relative' }}>
                    <Image
                      src="/logo-nav.png"
                      alt="SBC Logo"
                      fill
                      className="object-contain"
                    />
                  </div>
                  <div>
                    <h1 className="text-lg font-bold">SBC</h1>
                    <p className="text-xs text-muted-foreground">上海生物信息中心</p>
                  </div>
                </div>
                <Button variant="ghost" size="icon" onClick={() => setIsOpen(false)}>
                  <X className="h-5 w-5" />
                </Button>
              </div>
              <div className="mb-6">
                <SimpleGlobalSearch />
              </div>
              <nav className="space-y-4">
                {navigation.map((item) => {
                  const Icon = item.icon
                  return (
                    <Link
                      key={item.name}
                      href={item.href}
                      onClick={() => setIsOpen(false)}
                      className="flex items-center gap-3 text-lg font-medium transition-colors hover:text-primary"
                    >
                      {Icon && <Icon className="h-5 w-5" />}
                      {item.name}
                    </Link>
                  )
                })}
              </nav>
            </SheetContent>
          </Sheet>
        </div>
      </div>
    </header>
  )
}
